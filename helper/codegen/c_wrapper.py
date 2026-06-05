from .api.utils import CodeGenerator, Serializer, Indentable
from .api.fortran import Intent

INDENT = 4
TYPE_CONVERSION_SUFFIX = "_f"


class C_Wrapper_Serializer(Serializer):
    def Fortran_Type(self, spec):
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
                    case _:
                        raise ValueError(f"Unsupported conversion for '{self.name}'")

                if len(self.dimension) > 0 or self.name == "character":
                    string += ", allocatable"
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

        if len(self.dimension) > 0:
            string += format(self.dimension, spec)

        return Indentable(string)

    def DocList(self, spec):
        match self.type:
            case "argument":
                formatted = f"!!{"\n!!".join(self.doc_list)}"
            case "procedure" | "module":
                formatted = f"!>{"\n!|".join(self.doc_list)}"
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
        pass

    def Procedure_Arguments(self, spec):
        match spec:
            case "arglist":
                return Indentable(", ".join(arg.name + arg.type.needs_conversion * TYPE_CONVERSION_SUFFIX for arg in self))
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

    def Procedure(self, spec):
        match spec:
            case "call":
                call = f"{self.name}({format(self.args, "arglist")})"
                if self.retvar is None:
                    call = "call " + call
                else:
                    call = f"{self.retvar} = {call}"
                return Indentable(call)
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

    def C_Wrapper_Argument(self, spec):
        type_conversion_name = self.name + TYPE_CONVERSION_SUFFIX
        match spec:
            case "dummy":
                arg_str = f"{format(self.type, "dummy")}, target :: {self.name}"
                for line in self.doc:
                    arg_str += f"\n    !!{line}"
            case "locals_type_conversion":
                arg_str = f"{format(self.type, "locals_type_conversion")} :: {type_conversion_name}"
            case "type_conversion_inputs":
                if self.type.intent is Intent.OUT:
                    if self.type.name == "character":
                        arg_str = f"M_ALLOCATE(character(len={self.type.len}) :: {type_conversion_name}{format(self.type.dimension, "tuple")})"
                    else:
                        arg_str = f"M_ALLOCATE({type_conversion_name}{format(self.type.dimension, "tuple")})"
                else:
                    match self.type.name:
                        case "logical":
                            arg_str = f"call logical_as_c_int({self.name}, {type_conversion_name})"
                        case "character":
                            ndims = len(self.type.dimension)
                            if ndims == 0:
                                arg_str = f"M_ALLOCATE(character(len={self.type.len}) :: {type_conversion_name}{format(self.type.dimension, "tuple")})"
                                arg_str += f"\ncall char_as_c_char({self.name}, {type_conversion_name})"
                            elif ndims < 3:
                                arg_str = f"call string_as_c_char_{ndims}d({self.name}, {type_conversion_name})"
                        case _:
                            raise ValueError(f"Unsupported conversion for '{self.name}'")
            case "type_conversion_outputs":
                if self.type.intent is not Intent.IN:
                    match self.type.name:
                        case "logical":
                            arg_str = f"call c_int_as_logical({type_conversion_name}, {self.name})"
                        case "character":
                            ndims = len(self.type.dimension)
                            if ndims == 0:
                                arg_str = f"call c_char_as_char({self.name}, {type_conversion_name})"
                            elif ndims < 3:
                                arg_str = f"call c_char_{ndims}d_as_string({self.name}, {type_conversion_name})"
                        case _:
                            raise ValueError(f"Unsupported conversion for '{self.name}'")
                else:
                    raise ValueError(f"'{self.name}' is not output")
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")
        return Indentable(arg_str)

    def C_Wrapper_Arguments(self, spec):
        match spec:
            case "dummy":
                formatted_dim_args = "\n".join(format(arg, "dummy") for arg in self if arg.is_dim_arg)
                formatted_other = "\n".join(format(arg, "dummy") for arg in self if not arg.is_dim_arg)
                formatted = formatted_dim_args + "\n" + formatted_other
            case "locals_type_conversion" | "type_conversion_inputs":
                formatted = "\n".join(format(arg, spec) for arg in self if arg.type.needs_conversion)
            case "type_conversion_outputs":
                formatted = "\n".join(format(arg, spec) for arg in self if arg.type.needs_conversion and arg.type.intent is not Intent.IN)
            case "arglist":
                formatted = ", ".join(arg.name for arg in self)
            case "null_validation":
                formatted = "M_CHECK_IERR_NON_NULL"
                for arg in self:
                    if arg.name != "ierr":
                        formatted += f"\nM_CHECK_NON_NULL({arg.name})"
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

        return Indentable(formatted)

    def C_Wrapper(self, spec):
        orig_proc = self.orig_procedure
        module_name = orig_proc.parent.name
        ford_link = f"[[{module_name}(module):{orig_proc.name}({orig_proc.type})]]"
        doc = [f" C-wrapper for {ford_link}"] + self.doc

        wrapper = f"""{format(doc, "subroutine")}
subroutine {self.name}({format(self.arguments, "arglist")}) bind(C, name="{self.name}")
    use {module_name}, only: {self.orig_procedure.name}
{format(self.arguments, "dummy") >> INDENT}
{format(self.arguments, "locals_type_conversion") >> INDENT}
{format(self.arguments, "null_validation") >> INDENT}
{format(self.arguments, "type_conversion_inputs") >> INDENT}
{format(orig_proc, "call") >> INDENT}
{format(self.arguments, "type_conversion_outputs") >> INDENT}
end subroutine {self.name}"""

        return Indentable(wrapper)

    def C_Wrapper_Module(self, spec):
        doc = [f" Module for C-wrappers for [[{self.orig_module.name}(module)]]"] + self.doc
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

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

{"\n\n".join(format(c_wrapper) >> INDENT for c_wrapper in self)}

end module {self.name}
#endif"""
