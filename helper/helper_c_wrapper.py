import re

# name of the core subroutine
core_name = "sort_character"

# name of the module that contains the core subroutine
module_name = "tox_sorting"

# <--- Core subroutine arguments as a string --->
core_args = """
integer, intent(in) :: n, strlen
character(len=strlen), intent(in) :: array(n)
integer, intent(inout) :: perm(n), stack_left(n), stack_right(n)

"""

# <--- Define correct order of arguments that the subroutine receives --->
arg_order = [
    "array", "perm", "stack_left", "stack_right"
]

#-------------------------------------------------------------------------------------

type_map = {
    "integer": "integer(c_int)",
    "real(real64)": "real(c_double)",
    "real(8)": "real(c_double)"
}

def process_line(line):
    line = line.split('!')[0].strip()
    if not line:
        return []
    m = re.match(r'(\w+(?:\([^)]+\))?),?\s*(intent\([^)]+\))?\s*::\s*(.+)', line)
    if not m:
        return []
    ftype, intent, vars_str = m.groups()
    ftype_c = type_map.get(ftype.strip(), ftype.strip())
    intent_str = f", {intent}" if intent else ""
    decls = []
    for var in vars_str.split(','):
        var = var.strip()
        name = var.split('(')[0].strip()
        # Arrays: use (*) for C
        if '(' in var:
            decl = f"{ftype_c}{intent_str}, target :: {name}(*)"
        else:
            decl = f"{ftype_c}{intent_str}, value :: {name}"
        decls.append((name, decl))
    return decls

# Parse all declarations into a dict
decl_dict = {}
for line in core_args.strip().splitlines():
    for name, decl in process_line(line):
        decl_dict[name] = decl

# Build argument list and declarations in the desired order
subr_args = ", ".join(arg_order)
call_args = ", ".join(arg_order)
decl_lines = [decl_dict[name] for name in arg_order]

# Output wrapper
wrapper = []
wrapper.append(f"subroutine {core_name}_c({subr_args}) bind(C, name=\"{core_name}_c\")")
wrapper.append("  use iso_c_binding")
wrapper.append(f"  use {module_name}")
for decl in decl_lines:
    wrapper.append("  " + decl)
wrapper.append("")
wrapper.append(f"  call {core_name}({call_args})")
wrapper.append("end subroutine " + core_name + "_c")

print("\n".join(wrapper))