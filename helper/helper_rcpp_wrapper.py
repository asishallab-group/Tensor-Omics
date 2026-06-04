import re
import argparse
from textwrap import indent


def _preprocess_fortran(code):
    """Remove comments (but keep lines with special markers) and join continuations.

    This normalizes ampersand continuations so the parser can work on single logical lines.
    """
    lines = code.splitlines()
    out = []
    cont = ''
    for raw in lines:
        line = raw.rstrip()
        # Remove inline comment after '!' unless it's part of a marked comment starting with '!|'
        if '!|' in line:
            # keep the comment starting at '!|' for possible hints, but strip other trailing comments
            idx = line.find('!|')
            code_part = line[:idx]
            comment = line[idx:]
        else:
            # drop any inline comment
            idx = line.find('!')
            if idx != -1:
                code_part = line[:idx]
            else:
                code_part = line
            comment = ''

        code_part = code_part.rstrip()
        # Handle free-form Fortran continuation character '&' at end
        if code_part.endswith('&'):
            cont += code_part[:-1].rstrip() + ' '
        else:
            full = (cont + code_part).strip()
            cont = ''
            if full:
                out.append(full + (' ' + comment if comment else ''))
    return '\n'.join(out)


def generar_wrappers_fortran_a_c_rcpp(codigo_fortran, nombre_rcpp_func_prefix="tox_", nombre_rcpp_func_suffix="_rcpp"):
    """Parse a Fortran subroutine with bind(C) and generate a C header + Rcpp wrapper.

    The function attempts to be generic: it recognizes iso_c_binding types, intent(in/out/inout),
    arrays with dimensions, 'value' scalars, and maps common names for size variables to Rcpp expressions.

    Returns: (header_c_str, wrapper_rcpp_str)
    """
    code = _preprocess_fortran(codigo_fortran)

    # Extract subroutine name and bind(C name) (name may be optional)
    sub_pat = re.compile(r"subroutine\s+(\w+)\s*\((.*?)\)\s*bind\(C(?:,\s*name\s*=\s*\"(.*?)\")?\)", re.IGNORECASE | re.DOTALL)
    m = sub_pat.search(code)
    if not m:
        return "Error: No se encontró la declaración de subrutina 'bind(C)'.", ""

    subroutine_name = m.group(1)
    arglist = m.group(2)
    bind_name = m.group(3) or subroutine_name

    # Split arguments considering parentheses
    argumentos_brutos = re.split(r",\s*(?![^()]*\))", arglist)

    # Parse declarations line-by-line to avoid accidental multi-line captures
    decls = []
    for line in code.splitlines():
        s = line.strip()
        if re.match(r'^(integer|real|logical|character)\b', s, re.IGNORECASE):
            parts = s.split('::', 1)
            if len(parts) == 2:
                left = parts[0].strip()
                rights = parts[1].strip()
                decls.append((left, rights))

    # Build a map name->declaration
    decl_map = {}
    for left, rights in decls:
        left = left.strip()
        rights = rights.strip()
        for name in [n.strip() for n in re.split(r",\s*(?![^()]*\))", rights)]:
            arrm = re.match(r"(\w+)\s*\((.*?)\)", name)
            if arrm:
                nm = arrm.group(1)
                dims = arrm.group(2)
                decl_map[nm] = (left, nm, dims)
            else:
                decl_map[name] = (left, name, None)

    # Type mappings: Fortran iso C binding -> C type, and Rcpp containers
    tipo_c_map = {
        'integer(c_int)': 'int',
        'integer(c_long)': 'long',
        'real(c_double)': 'double',
        'logical(c_bool)': 'int',  # use int for C header compatibility
        'character(c_char)': 'char',
    }

    rcpp_type_map = {
        'integer': 'int',
        'real': 'double',
    }

    # Helpers
    def normalize_type(t):
        return re.sub(r"\s+", "", t.lower())

    args_c = []
    rcpp_inputs = []
    rcpp_body_vars = []
    rcpp_call_args = []
    rcpp_return_vars = []
    c_call_args = []
    # Body lines collected during generation (declare early to allow filename handling)
    body_lines = []

    # For size inference heuristics
    size_inference = {
        'n_axes': "expression_vectors.nrow()",
        'n_vectors': "expression_vectors.ncol()",
        'n_selected_vectors': "sum(vector_selection)",
        'n_selected_axes': "sum(axis_selection)",
    }

    # No hardcoded explicit inputs: we'll infer Rcpp types and (optionally) friendly names

    def friendly_name(fortran_name):
        n = fortran_name.lower()
        # heuristics to produce clearer Rcpp arg names
        if 'exp_vec' in n or ('vec' in n and 'select' in n) or ('selection' in n and 'vec' in n) or n.endswith('_index'):
            return 'vector_selection'
        if 'axis' in n and 'select' in n:
            return 'axis_selection'
        return fortran_name

    # Map fortran arg order -> we will iterate in the order given in arglist
    # Pre-check for filename_ascii pattern
    argumentos_names = [a.strip() for a in argumentos_brutos]
    has_filename_ascii = 'filename_ascii' in argumentos_names and 'fn_len' in argumentos_names
    # body marker to avoid duplicate injected lines
    body_marker_set = set()

    for arg in argumentos_brutos:
        name = arg.strip()
        if not name:
            continue
        # Detect filename_ascii pattern: if present we'll accept a single `std::string filename` in R
        # and construct `filename_ascii` + `fn_len` in the wrapper body.
        # We don't add `filename_ascii` as an explicit R parameter.
        # Record presence for later use.
        # (This is a simple heuristic; keeps generator conservative.)
        #
        # Note: argumentos_brutos is the Fortran arglist in order.
        # If a Fortran subroutine uses both `filename_ascii` and `fn_len`, we'll synthesize `filename`.
        pass
        if name not in decl_map:
            # if missing declaration, assume integer input scalar
            decl = ('integer(c_int), intent(in), value', name, None)
        else:
            decl = decl_map[name]

        left, nm, dims = decl
        left_parts = [p.strip() for p in left.split(',') if p.strip()]
        ftype = left_parts[0]
        intent = None
        is_value = False
        is_target = False
        for p in left_parts[1:]:
            if p.lower().startswith('intent'):
                intent = re.search(r"intent\((in|out|inout)\)", p, re.IGNORECASE).group(1).lower() if re.search(r"intent\((in|out|inout)\)", p, re.IGNORECASE) else None
            if p.lower() == 'value':
                is_value = True
            if p.lower() == 'target':
                is_target = True

        norm_ft = normalize_type(ftype)
        ctype = tipo_c_map.get(norm_ft, ftype)

        is_scalar = (dims is None)

        # Special handling: if this Fortran subroutine uses filename_ascii + fn_len,
        # expose a single `std::string filename` parameter to the R wrapper and
        # synthesize `filename_ascii` + `fn_len` inside the wrapper body.
        if has_filename_ascii and nm == 'filename_ascii':
            # Add R-level `std::string filename` if not already present
            if ('filename', 'std::string') not in rcpp_inputs:
                rcpp_inputs.append(('filename', 'std::string'))
            # Inject filename conversion body lines once
            if 'filename_ascii_build' not in body_marker_set:
                body_lines.append('int fn_len = (int)filename.size();')
                body_lines.append('IntegerVector filename_ascii(fn_len);')
                body_lines.append('for (int i = 0; i < fn_len; ++i) filename_ascii[i] = (int)filename[i];')
                body_marker_set.add('filename_ascii_build')
            # Pass the constructed buffer and length to the C call
            c_call_args.append('filename_ascii.begin()')
            c_call_args.append('fn_len')
            # Still include the underlying C prototype (args_c) below as normal

        # Build C arg
        if is_scalar:
            if intent == 'in' and is_value:
                args_c.append(f"{ctype} {nm}")
            elif intent == 'in' and not is_value:
                # Fortran may pass scalars without value - treat as pointer
                args_c.append(f"{ctype}* {nm}")
            else:
                args_c.append(f"{ctype}* {nm}")
        else:
            # arrays passed as pointers
            args_c.append(f"{ctype}* {nm}")

        # Infer Rcpp input/output and C-call arguments (preserve Fortran order)
        if intent == 'in':
            if not is_scalar:
                # Two-dimensional arrays -> NumericMatrix (common pattern)
                if dims and ',' in dims:
                    r_name = friendly_name(nm)
                    r_type = 'NumericMatrix'
                    if (r_name, r_type) not in rcpp_inputs:
                        rcpp_inputs.append((r_name, r_type))
                    c_call_args.append(f"{r_name}.begin()")
                else:
                    # One-dimensional array -> Vector (Integer/Numeric)
                    r_type = 'IntegerVector' if 'integer' in norm_ft or 'logical' in norm_ft else 'NumericVector'
                    r_name = friendly_name(nm)
                    if (r_name, r_type) not in rcpp_inputs:
                        rcpp_inputs.append((r_name, r_type))
                    c_call_args.append(f"{r_name}.begin()")
            else:
                # scalar input: if it's a size we will compute it in Rcpp (don't make it a parameter)
                if nm in size_inference:
                    # body var will be added later; pass the local var name to C call
                    c_call_args.append(nm)
                else:
                    r_type = 'int' if 'integer' in norm_ft else 'double'
                    # include simple scalars as parameters
                    if (nm, r_type) not in rcpp_inputs:
                        rcpp_inputs.append((nm, r_type))
                    if is_value:
                        c_call_args.append(nm)
                    else:
                        c_call_args.append(f"&{nm}")

        # If the Fortran signature included fn_len we already handled passing it
        if has_filename_ascii and nm == 'fn_len':
            # skip adding fn_len as a separate R parameter or c_call arg (already appended)
            continue

        elif intent in ('out', 'inout'):
            if is_scalar:
                rcpp_return_vars.append({'name': nm, 'is_scalar': True, 'ctype': ctype})
                c_call_args.append(f"&{nm}")
            else:
                size_token = dims.split(',')[0].strip() if dims else '0'
                size_expr = size_inference.get(size_token, size_token)
                rcpp_type = 'NumericVector' if 'real' in norm_ft else 'IntegerVector'
                rcpp_return_vars.append({'name': nm, 'is_scalar': False, 'type': rcpp_type, 'size': size_expr})
                c_call_args.append(f"{nm}.begin()")

        else:
            # Unknown intent -> treat as input
            if not is_scalar:
                r_type = 'NumericVector' if 'real' in norm_ft else 'IntegerVector'
                r_name = friendly_name(nm)
                if (r_name, r_type) not in rcpp_inputs:
                    rcpp_inputs.append((r_name, r_type))
                c_call_args.append(f"{r_name}.begin()")
            else:
                r_type = 'int' if 'integer' in norm_ft else 'double'
                if (nm, r_type) not in rcpp_inputs:
                    rcpp_inputs.append((nm, r_type))
                c_call_args.append(nm)

    # Create header C (multi-line for readability)
    if args_c:
        header_c = f"void {bind_name}(\n    " + ",\n    ".join(args_c) + "\n);"
    else:
        header_c = f"void {bind_name}(void);"

    # Build Rcpp function signature: include only explicit inputs in a convenient order
    # Prefer common order: expression_vectors, vector_selection, axis_selection, others
    preferred_order = ['expression_vectors', 'exp_vecs_selection_index', 'exp_vecs_selection', 'vector_selection', 'axes_selection', 'axis_selection']
    inputs_sorted = []
    seen = set()
    for pname in preferred_order:
        for (n, t) in rcpp_inputs:
            if n == pname and n not in seen:
                inputs_sorted.append((n, t))
                seen.add(n)
    for (n, t) in rcpp_inputs:
        if n not in seen:
            inputs_sorted.append((n, t))
            seen.add(n)

    # Build function parameter string
    params_decl = ',\n'.join([f"                                      {t} {n}" for n, t in inputs_sorted])
    if not params_decl:
        params_decl = ''

    # Body: declare size vars and outputs
    # Size vars from heuristics
    for _, t in inputs_sorted:
        pass

    body_lines = []
    # Declare inferred size variables if they are used in rcpp_return_vars
    # Also declare any scalar input parameters that were added as plain types (int/double) - no extra action needed
    # For example n_axes/n_vectors/n_selected_vectors
    for arg in argumentos_brutos:
        name = arg.strip()
        if name in size_inference:
            body_lines.append(f"int {name} = {size_inference[name]};")

    # Declare outputs
    for var in rcpp_return_vars:
        if var.get('is_scalar'):
            body_lines.append(f"{var['ctype']} {var['name']} = 0;")
        else:
            body_lines.append(f"{var['type']} {var['name']}({var['size']});")

    # For any inout scalars that were represented as pointer args, ensure local declarations exist
    # (already handled above as return scalars)

    # Build the call line
    call_args = ', '.join(c_call_args)

    # Decide return statement. We will include:
    #  - array outputs (e.g., tissue_versatilities)
    #  - inferred size vars (n_selected_vectors, n_selected_axes) when computed
    #  - scalar outputs (e.g., ierr)
    used_size_vars = [n for n in size_inference.keys() if n in [a.strip() for a in argumentos_brutos]]

    # Only return Fortran outputs (intent(out) or inout). rcpp_return_vars already contains these.
    if len(rcpp_return_vars) == 0:
        return_type = 'SEXP'
        return_stmt = '    return R_NilValue;'
    elif len(rcpp_return_vars) == 1:
        rv = rcpp_return_vars[0]
        # Single output: return its Rcpp type directly (NumericVector/IntegerVector) or scalar C type (int/double)
        if rv.get('is_scalar'):
            return_type = rv['ctype']
        else:
            return_type = rv['type']
        return_stmt = f"    return {rv['name']};"
    else:
        return_type = 'List'
        items = []
        # arrays and scalars together
        for var in rcpp_return_vars:
            items.append(f"        Named(\"{var['name']}\") = {var['name']}")

        return_stmt = "    return List::create(\n" + ",\n".join(items) + "\n    );"

    # Normalize bind name for wrapper (strip common prefixes/suffixes)
    def normalize_bind_name(name):
        n = name
        # strip trailing _c, _f, _f90
        n = re.sub(r"(_c|_f90|_f)$", "", n, flags=re.IGNORECASE)
        # strip common leading verbs
        n = re.sub(r"^(compute_|calculate_|calc_|get_|set_)", "", n, flags=re.IGNORECASE)
        return n

    nombre_rcpp = nombre_rcpp_func_prefix + normalize_bind_name(bind_name) + nombre_rcpp_func_suffix

    params_sig = ',\n'.join([f"{t} {n}" for n, t in inputs_sorted])
    if params_sig:
        params_sig = params_sig

    body = []
    if body_lines:
        body += [l for l in body_lines]

    body.append(f"{bind_name}({call_args});")
    if return_stmt:
        body.append(return_stmt)

    body_str = '\n    '.join(body)

    # Build the final wrapper text
    # Build params header: create signature with defaults if none
    if inputs_sorted:
        params_text = ',\n'.join([f"{t} {n}" for n, t in inputs_sorted])
        wrapper = f"// [[Rcpp::export]]\n{return_type} {nombre_rcpp}({params_text}) {{\n    {body_str}\n}}"
    else:
        wrapper = f"// [[Rcpp::export]]\n{return_type} {nombre_rcpp}() {{\n    {body_str}\n}}"

    # Pretty indenting
    wrapper = re.sub(r"\n\s+", "\n    ", wrapper)

    return header_c, wrapper


def generate_wrappers_fortran_to_c_rcpp(fortran_code, rcpp_prefix="tox_", rcpp_suffix="_rcpp"):
    """English wrapper: call the internal Fortran->C/Rcpp generator.

    This is a thin English-named API that delegates to the parser implementation
    (keeps backward compatibility with the original Spanish-named function).

    Args:
        fortran_code (str): Fortran subroutine source.
        rcpp_prefix (str): Prefix for the generated Rcpp function name.
        rcpp_suffix (str): Suffix for the generated Rcpp function name.

    Returns:
        tuple: (c_header_str, rcpp_wrapper_str)
    """
    # Delegate to existing implementation (keeps parsing logic unchanged)
    return generar_wrappers_fortran_a_c_rcpp(fortran_code, rcpp_prefix, rcpp_suffix)


def _demo_examples():
    # Two examples: the one provided and a small scalar example
    fortran_example_code = """
subroutine distance_to_centroid_c(n_genes, n_families, genes, centroids, & 
                                  gene_to_fam, distances, d) bind(C, name="distance_to_centroid_c")
  use iso_c_binding, only : c_int, c_double
  use tox_euclidean_distance
  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes
  !| Total number of gene families
  integer(c_int), intent(in), value :: n_families
  !| Expression vector dimension
  integer(c_int), intent(in), value :: d
  !| Gene expression matrix (d × n_genes), column-major
  real(c_double), intent(in), target :: genes(d, n_genes)
  !| Family centroid matrix (d × n_families), column-major
  real(c_double), intent(in), target :: centroids(d, n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  !| Output distances array
  real(c_double), intent(out), target :: distances(n_genes)
  
  call distance_to_centroid(n_genes, n_families, genes, centroids, &
                            gene_to_fam, distances, d)
end subroutine distance_to_centroid_c

"""

    header_c, wrapper_rcpp = generate_wrappers_fortran_to_c_rcpp(
        fortran_example_code,
        rcpp_prefix="tox_",
        rcpp_suffix="_rcpp"
    )

    print("## C header generated")
    print(header_c)
    print('\n## Rcpp (C++) wrapper generated')
    print(wrapper_rcpp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate C header and Rcpp wrapper from a Fortran bind(C) subroutine')
    parser.add_argument('file', nargs='?', help='Fortran file (optional). If omitted, demo examples are shown.')
    args = parser.parse_args()
    if not args.file: 
        _demo_examples()
    else:
        with open(args.file, 'r') as f:
            code = f.read()
        header, wrapper = generate_wrappers_fortran_to_c_rcpp(code)
        print(header)
        print()
        print(wrapper)