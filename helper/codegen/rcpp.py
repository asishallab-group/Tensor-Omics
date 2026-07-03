from .api.utils import CodeGenerator, Serializer, Indentable
from .api.c_wrapper import C_Wrapper_Modules
from .api.fortran import Intent, Dimension
import re


INDENT = 4
DIM_ARG_RE = re.compile(r"^[^\d\*]")


class Rcpp_Serializer(Serializer):
    def Fortran_Type(self, spec):
        match spec:
            case "ctype":
                match self.name:
                    case "integer":
                        string = "int"
                    case "logical":
                        string = "int"
                    case "real":
                        string = "double"
                    case "complex":
                        string = "double _Complex"
                    case "character":
                        string = "char"
            case "rcpptype":
                ndim = self.ndim
                if ndim == 0 and self.name != "logical":
                    return format(self, "ctype")
                else:
                    match self.name:
                        case "integer":
                            string = "Integer"
                        case "logical":
                            string = "Logical"
                            ndim = 3
                        case "real":
                            string = "Numeric"
                        case "complex":
                            string = "Complex"
                        case "character":
                            ndim -= 1
                            string = "String"
                    if ndim == 1 or ndim > 2:
                        string += "Vector"
                    elif ndim == 2:
                        string += "Matrix"
            case "arglist":
                string = f"{format(self, "rcpptype")}"
                if self.ndim > 0:
                    string += "&"
                if self.intent is Intent.IN:
                    string = f"const {string}"
            case "extern":
                ndim = self.ndim
                string = f"{format(self, "ctype")}*"
                if self.intent is Intent.IN:
                    string = f"const {string}"
            case "init_val":
                match self.name:
                    case "integer":
                        string = "0"
                    case "logical":
                        string = "0"
                    case "real":
                        string = "0.0"
                    case "complex":
                        string = "0.0"
                    case "character":
                        string = '"\\0"'

        return Indentable(string)

    def DocList(self, spec):
        raise NotImplementedError()

    def Dimension(self, spec):
        match spec:
            case "flat_size":
                string = " * ".join(self)

        return Indentable(string)

    def Intent(self, spec):
        raise NotImplementedError()

    def Procedure(self, spec):
        raise NotImplementedError()

    def Procedure_Argument(self, spec):
        raise NotImplementedError()

    def Procedure_Arguments(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Argument(self, spec):
        error_handling = self.parent.orig_procedure.parent.error_handling
        ERR_NAN_INF = error_handling.error_codes["ERR_NAN_INF"][0]
        def return_error(code): return f'return List::create(Named("ierr") = {code});'

        string = ""
        match spec:
            case "extern":
                string = f"{format(self.type, spec)} {self.name}"
            case "arglist":
                if self.type.intent is not Intent.OUT and len(self.is_dim_arg_for) == 0:
                    default_value = self.default_value or ""
                    if self.optional:
                        default_value = f" = {default_value}"
                    string = f"{format(self.type, spec)} {self.name}{default_value}"
            case "return":
                if self.type.intent is Intent.OUT and not self.is_temporary:
                    string = f'Named("{self.name}") = {self.name}'
            case "call":
                # scalars as pointer
                if self.type.ndim == 0:
                    string = f"&{self.name}"
                # work arrays are std:vector objects
                elif self.is_temporary:
                    string = f"{self.name}.data()"
                # type conversions always create the pointer during conversion
                elif self.type.needs_conversion:
                    string = f"{self.name}_p"
                # non-temporary arrays are Rcpp objects
                else:
                    string = f"{self.name}.begin()"
            case "allocate":
                if self.type.intent is Intent.OUT:
                    if self.type.ndim == 0:
                        string = f"{format(self.type, "ctype")} {self.name} = {format(self.type, "init_val")};"
                    else:
                        # for temporaries, create std::vector, but also for character, as it needs to be converted
                        if self.is_temporary or self.type.name == "character":
                            if self.type.name == "character":
                                init_val = format(self.type, "init_val")
                                init_val = "" if not init_val else f", {init_val}"
                            else:
                                init_val = ""
                            string = f'std::vector<{format(self.type, "ctype")}> {self.name}({format(self.type.dimension, "flat_size")}{init_val});'
                        else:
                            string = f'{format(self.type, "rcpptype")} {self.name}({format(self.type.dimension, "flat_size")});'
            case "type_conversion_inputs":
                if self.type.needs_conversion and self.type.intent is not Intent.OUT:
                    pointer = self.name + "_p"
                    if self.type.name == "logical":
                        string = f"""{format(self.type, "extern")} {pointer} = {self.name}.begin();
// Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
for (int i = 0; i < {self.name}.size(); ++i) {{
    if ({pointer}[i] == NA_LOGICAL) {{
        {return_error(ERR_NAN_INF)}
    }}
}}"""
                    else:
                        if self.type.ndim == 1:
                            str_len = self.type.dimension[0]
                            str_len_init = ""
                            if str_len < "0" or str_len > "9":
                                str_len_init = f"\nint {str_len} = std::strlen({pointer});"

                            string = f"""if ({self.name} == NA_STRING) {{
    {return_error(ERR_NAN_INF)}
}}
{format(self.type, "extern")} {pointer} = {self.name}.get_cstring();{str_len_init}"""
                        else:
                            sexp_pointer = self.name + "_SEXP_p"
                            size = f"{self.name}.size()" if self.type.ndim != 3 else f"{self.name}.ncol() * {self.name}.nrow()"
                            max_str_len = self.type.dimension[0]
                            max_str_len_init = ""
                            if max_str_len < "0" or max_str_len > "9":
                                max_str_len_init = f"""
int {max_str_len} = 0;
for (int i = 0; i < {size}; ++i)
    {max_str_len} = std::max({max_str_len}, (int)std::strlen(CHAR({pointer}[i])));
"""
                            string = f"""SEXP* {sexp_pointer} = {self.name}.begin();{max_str_len_init}
std::vector<{format(self.type, "ctype")}> {self.name}_c({size});
{format(self.type, "extern")} {pointer} = {self.name}_c.data();
// Check for NA values and return ERR_NAN_INF code, and convert
for (int i = 0; i < {size}; ++i) {{
    if ({sexp_pointer}[i] == NA_STRING) {{
        {return_error(ERR_NAN_INF)}
    }} else {{
        const char* str = CHAR({sexp_pointer}[i]);
        int len = std::min((int)std::strlen(str), {max_str_len});
        std::memcpy({pointer} + i * {max_str_len}, str, len);
    }}
}}"""
            case "type_conversion_outputs":
                string = ""

        return Indentable(string)

    def C_Wrapper_Arguments(self, spec):
        match spec:
            case "extern" | "arglist" | "return" | "call":
                string = ",\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "allocate" | "type_conversion_inputs":
                string = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "dim_args":
                seen_dim_args = set()
                string = ""
                for arg in self:
                    if arg.type.ndim > 0:
                        dim_args = []
                        dimension = arg.type.dimension if arg.type.name != "character" else arg.type.dimension[1:]
                        for i, dim_arg_name in enumerate(dimension):
                            if DIM_ARG_RE.match(dim_arg_name) is not None and dim_arg_name not in seen_dim_args:
                                seen_dim_args.add(dim_arg_name)
                                dim_args.append((i, dim_arg_name))

                        if dim_args:
                            name = f"{arg.name}_shape"
                            string += f"""
IntegerVector {name} = {arg.name}.attr("dim");
{"\n".join(f"{dim_arg} = {name}[{i}]" for i, dim_arg in dim_args)}"""

        return Indentable(string)

    def C_Wrapper(self, spec):
        if spec == "extern":
            formatted = f"""void {self.name}(
{format(self.arguments, spec) >> INDENT}
);"""
        else:
            formatted = f"""List {self.stripped_name + "_rcpp"}(
{format(self.arguments, "arglist") >> INDENT}
) {{

{format(self.arguments, "dim_args") >> INDENT}

{format(self.arguments, "type_conversion_inputs") >> INDENT}

{format(self.arguments, "allocate") >> INDENT}

    {self.name}(
{format(self.arguments, "call") >> INDENT * 2}
    );

//{{format(self.arguments, "type_conversion_outputs")}}

    return List::create(
{format(self.arguments, "return") >> 2 * INDENT}
    );
}}"""
        return Indentable(formatted)

    def C_Wrapper_Module(self, spec):
        formatted = "\n\n".join(format(c_wrapper, spec) for c_wrapper in self.c_wrappers)
        match spec:
            case "extern":
                formatted = f"""extern "C" {{
{Indentable(formatted) >> INDENT}
}}"""

        return formatted

    def C_Wrapper_Modules(self, spec):
        return "\n\n".join(format(mod) for mod in self)

    @classmethod
    def dump(cls, c_wrapper_modules: C_Wrapper_Modules, out_dir: str = "rcpp/tensor_omics"):
        from pathlib import Path
        from os import mkdir

        c_wrapper_modules.use(cls)
        out_dir = Path(out_dir)
        if not out_dir.parent.is_dir():
            raise ValueError(f"out_dir's directory '{out_dir.parent}' does not exist")

        if out_dir.is_dir():
            from shutil import rmtree
            rmtree(out_dir)

        mkdir(out_dir)

        for module in c_wrapper_modules:
            with open(out_dir.joinpath(module.stripped_name + ".cpp"), "w") as module_file:
                module_file.write(f"""#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

{format(module, "extern")}

{format(module)}""")
