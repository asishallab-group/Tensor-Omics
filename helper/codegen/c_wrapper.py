from .api.utils import CodeGenerator, Serializer, Indentable
from .api.fortran import Intent
from .api.c_wrapper import C_Wrapper_Modules

INDENT = 4
TYPE_CONVERSION_SUFFIX = "_f"
MAPPED_MODE_SUFFIX = "_int_f"


class C_Wrapper_Serializer(Serializer):
    def Fortran_Type(self, spec):
        dimension = self.dimension

        match spec:
            case "dummy":
                match self.name:
                    case "logical":
                        string = f"integer(c_int)"
                    case "character":
                        string = f"character(len=1, kind=c_char)"
                    case "integer":
                        string = f"integer(c_int)"
                    case "real":
                        string = f"real(c_double)"
                    case "complex":
                        string = f"real(c_double_complex)"
                    case _:
                        raise ValueError(f"'{self.name}' not supported for C formatting")

                string += format(self.intent)
            case "locals_type_conversion":
                match self.name:
                    case "logical":
                        string = f"logical"
                    case "character":
                        string = f"character(len=:)"
                        dimension = self.dimension[1:]
                    case _:
                        raise ValueError(f"Unsupported conversion for '{self.name}'")

                if len(self.dimension) > 0 or self.name == "character":
                    string += ", allocatable"
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

        if len(dimension) > 0:
            string += format(dimension, spec)

        return Indentable(string)

    def DocList(self, spec):
        match self.type:
            case "argument":
                formatted = f"!! {"\n!! ".join(self.doc_list)}"
            case "procedure" | "module":
                formatted = f"!> {"\n!| ".join(self.doc_list)}"
        return Indentable(formatted)

    def Dimension(self, spec):
        if len(self) == 0:
            return ""
        else:
            match spec:
                case "dummy":
                    return f", dimension({", ".join(self)})"
                case "locals_type_conversion":
                    return f", dimension({", ".join(":" for i in self)})"
                case "tuple":
                    return f"({", ".join(self)})"
                case _:
                    raise ValueError(f"Unsupported format spec '{spec}'")

    def Intent(self, spec):
        return f", intent({self.name.lower()})"

    def Procedure_Argument(self, spec):
        match spec:
            case "arglist":
                if self.is_mode_arg:
                    arg = f"{self.name} = {self.name + MAPPED_MODE_SUFFIX}"
                elif self.type.needs_conversion:
                    arg = f"{self.name} = {self.name + TYPE_CONVERSION_SUFFIX}"
                else:
                    arg = f"{self.name} = {self.name}"
                return Indentable(arg)
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

    def Procedure(self, spec):
        match spec:
            case "call":
                # arguments that will always be part of the call, either because non-optional or having a default value
                fixed_args = list(arg for arg in self.args if not arg.optional or arg.default_value is not None)

                # arguments whose presence depends on other args
                required_if_args = tuple(arg for arg in self.args if arg.doc_list.meta["required_if_mode"] is not None)
                required_if_args_targets = {}
                for required_if_arg in required_if_args:
                    meta = required_if_arg.doc_list.meta["required_if_mode"]
                    mode = meta.group("mode_var_name")
                    if mode not in required_if_args_targets:
                        required_if_args_targets[mode] = {}

                    mode_var = meta.group("mode_var")
                    if mode_var not in required_if_args_targets[mode]:
                        required_if_args_targets[mode][mode_var] = []

                    required_if_args_targets[mode][mode_var].append(required_if_arg)

                # optional outputs
                optional_outputs = tuple(arg for arg in self.args if arg.doc_list.meta["optional_output"] is not None)

                call_prefix = "call" if self.retvar is None else f"{self.retvar} ="

                def optional_variants(args, n_remaining):
                    if n_remaining > 0:
                        if len(args) == n_remaining:
                            call_without_optionals = "else\n" + (optional_variants([], 0) >> INDENT) + "\n"
                        current = args[0]
                        if type(current) is str:
                            mode = current
                            items = list(required_if_args_targets[mode].items())
                            variant = ""
                            while len(items) > 1:
                                mode_var, mode_args = items.pop()
                                condition = f"{mode + MAPPED_MODE_SUFFIX} == {mode_var}"
                                variant += f"""if ({mode + MAPPED_MODE_SUFFIX} == {mode_var}) then
{optional_variants(args[1:] + mode_args, n_remaining - 1) >> INDENT}
else """
                            mode_var, mode_args = items.pop()
                            condition = f"{mode + MAPPED_MODE_SUFFIX} == {mode_var}"
                            variant += f"""if ({mode + MAPPED_MODE_SUFFIX} == {mode_var}) then
{optional_variants(args[1:] + mode_args, n_remaining - 1) >> INDENT}
{call_without_optionals}end if"""
                        elif meta["optional_output"] is not None:
                            variant = f"""
if (present{current}) then
{optional_variants(args[1:] + [current], n_remaining - 1) >> INDENT}
else if (.not. present({current}))
{optional_variants(args[1:], n_remaining - 1) >> INDENT}
{call_without_optionals}end if
"""
                    else:
                        variant = f"""{call_prefix} {self.name}(&
{",&\n".join(format(arg, "arglist") >> INDENT for arg in fixed_args + args)}&
)"""

                    return Indentable(variant)

                optional_variant_args = list(optional_outputs)
                optional_variant_args.extend(required_if_args_targets.keys())

                call = optional_variants(optional_variant_args, len(optional_variant_args))
                return Indentable(call)
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

    def C_Wrapper_Argument(self, spec):
        type_conversion_name = self.name + TYPE_CONVERSION_SUFFIX
        arg_str = ""
        match spec:
            case "dummy":
                arg_str = f"{format(self.type, "dummy")}, target :: {self.name}\n"
                arg_str += format(self.doc_list, "argument") >> INDENT
            case "locals_type_conversion":
                if self.type.needs_conversion:
                    arg_str = f"{format(self.type, "locals_type_conversion")} :: {type_conversion_name}"
                    if self.mode_vars is not None:
                        arg_str += f"\ninteger(int32) :: {self.name}{MAPPED_MODE_SUFFIX}"
            case "type_conversion_inputs":
                if self.type.needs_conversion:
                    shape_arg = self.shape_arg
                    if shape_arg is None:
                        if self.type.name == "character":
                            shape = self.type.dimension[1:]
                        else:
                            shape = self.type.dimension
                        c_name = self.name
                    else:
                        shape = shape_arg.type.dimension
                        if self.type.name == "character":
                            c_name = f"{self.name}(:, 1:size({shape_arg.name}, kind=int32))"
                        else:
                            c_name = f"{self.name}(1:size({shape_arg.name}, kind=int32))"

                    if self.type.intent is Intent.OUT and len(self.type.dimension) > 0:
                        if self.type.name == "character":
                            arg_str = f"M_ALLOCATE(character(len={self.type.dimension[0]}) :: {type_conversion_name}{format(shape, "tuple")})"
                        else:
                            arg_str = f"M_ALLOCATE({type_conversion_name}{format(shape, "tuple")})"
                    else:
                        match self.type.name:
                            case "logical":
                                if len(self.type.dimension) > 0:
                                    arg_str = f"M_ALLOCATE({type_conversion_name}{format(shape, "tuple")})\n"
                                else:
                                    arg_str = ""
                                arg_str += f"call c_int_as_logical({c_name}, {type_conversion_name})"
                            case "character":
                                ndims = len(self.type.dimension)
                                if ndims == 0:
                                    arg_str = f"call c_char_as_char({c_name}, {type_conversion_name})"
                                elif ndims < 3:
                                    arg_str = f"call c_char_{ndims}d_as_string({c_name}, {type_conversion_name}, ierr)\nif (is_err(ierr)) return"
                                else:
                                    raise ValueError(f"String conversion for {ndims}D arrays not supported yet in tox conversions")
                            case _:
                                raise ValueError(f"Unsupported conversion for '{self.name}'")

                    if self.mode_vars is not None:
                        arg_str += f"\n\nselect case ({self.name}_f)"
                        for mode_name, module_name, mode_str in self.mode_vars:
                            arg_str += f"""\n    case ("{mode_str}")
            {self.name}{MAPPED_MODE_SUFFIX} = {mode_name}"""
                        arg_str += "\n    case default"
                        arg_str += "\n        call set_err(ierr, ERR_INVALID_INPUT)"
                        arg_str += "\n        return"
                        arg_str += "\nend select"

                    elif self.doc_list.meta["required_if_mode"] is not None:
                        required_if_mode = self.doc_list.meta["required_if_mode"]
                        arg_str = f"""if ({required_if_mode.group("mode_var_name") + MAPPED_MODE_SUFFIX} == {required_if_mode.group("mode_var")}) then
    M_CHECK_NON_NULL({self.name})
{Indentable(arg_str) >> INDENT}
end if"""
            case "type_conversion_outputs":
                if self.type.needs_conversion and self.type.intent is not Intent.IN:
                    shape_arg = self.shape_arg
                    if shape_arg is None:
                        c_name = self.name
                    else:
                        if self.type.name == "character":
                            c_name = f"{self.name}(:, 1:size({shape_arg.name}, kind=int32))"
                        else:
                            c_name = f"{self.name}(1:size({shape_arg.name}, kind=int32))"
                    match self.type.name:
                        case "logical":
                            arg_str = f"call logical_as_c_int({type_conversion_name}, {c_name})"
                        case "character":
                            ndims = len(self.type.dimension)
                            arg_str = f"call string_as_c_char_{ndims}d({type_conversion_name}, {c_name})"
                        case _:
                            raise ValueError(f"Unsupported conversion for '{self.name}'")
                else:
                    arg_str = ""
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")
        return Indentable(arg_str)

    def C_Wrapper_Arguments(self, spec):
        match spec:
            case "dummy":
                formatted_dim_args = "\n".join(format(arg, "dummy") for arg in self if len(arg.is_dim_arg_for) > 0)
                formatted_other = "\n".join(format(arg, "dummy") for arg in self if len(arg.is_dim_arg_for) == 0)
                formatted = formatted_dim_args + "\n" + formatted_other
            case "locals_type_conversion" | "mode_var_conversion":
                formatted = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "type_conversion_inputs":
                mode_vars = "\n".join(filter(bool, (format(arg, spec) for arg in self if arg.mode_vars is not None)))
                non_mode_vars = "\n".join(filter(bool, (format(arg, spec) for arg in self if arg.mode_vars is None)))
                formatted = mode_vars + "\n" + non_mode_vars
            case "type_conversion_outputs":
                formatted = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "arglist":
                formatted = ",&\n".join(arg.name for arg in self)
            case "null_validation":
                formatted = "M_CHECK_IERR_NON_NULL"
                for arg in self:
                    if arg.name != "ierr" and not arg.optional:
                        formatted += f"\nM_CHECK_NON_NULL({arg.name})"

            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

        return Indentable(formatted)

    def C_Wrapper(self, spec):
        orig_proc = self.orig_procedure
        module_name = orig_proc.parent.name
        ford_link = f"[[{module_name}(module):{orig_proc.name}({orig_proc.type})]]"
        doc = [f"summary: C-wrapper for {ford_link}"] + self.doc_list

        wrapper = f"""{format(doc, "subroutine")}
subroutine {self.name}(&
{format(self.arguments, "arglist") >> 2 * INDENT}&
        ) bind(C, name="{self.name}")
    use {module_name}, only: {self.orig_procedure.name}
    use {module_name}
{format(self.arguments, "dummy") >> INDENT}
{format(self.arguments, "locals_type_conversion") >> INDENT}
{format(self.arguments, "null_validation") >> INDENT}
{format(self.arguments, "type_conversion_inputs") >> INDENT}
{format(orig_proc, "call") >> INDENT}
{format(self.arguments, "type_conversion_outputs") >> INDENT}
end subroutine {self.name}"""

        return Indentable(wrapper)

    def C_Wrapper_Module(self, spec):
        doc = [f"summary: Module for C-wrappers for [[{self.orig_module.name}(module)]]"] + self.doc_list
        return f"""#ifndef NO_C_INTERFACE
#include <src/macros.h>

{format(doc, "module")}
module {self.name}
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

{"\n\n".join(format(c_wrapper) >> INDENT for c_wrapper in self)}

end module {self.name}
#endif"""

    @classmethod
    def dump(cls, c_wrapper_modules: C_Wrapper_Modules, out_dir: str):
        from pathlib import Path
        from os import mkdir

        c_wrapper_modules.use(cls)
        out_dir = Path(out_dir)
        if not out_dir.is_dir():
            raise ValueError("out_dir muste be a valid directory path")

        c_wrapper_dir = out_dir.joinpath("c_interface")
        if c_wrapper_dir.is_dir():
            from shutil import rmtree
            rmtree(c_wrapper_dir)

        mkdir(c_wrapper_dir)

        for module in c_wrapper_modules:
            module_file_path = c_wrapper_dir.joinpath(module.name + ".F90")

            with open(module_file_path, "w") as module_file:
                module_file.write(format(module))
