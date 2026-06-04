from .api.utils import CodeGenerator, Serializer, Indentable

INDENT = 4


class C_Wrapper_Serializer(Serializer):
    def Fortran_Type(self, spec):
        match self.name:
            case "logical":
                return f"integer(c_int)"
            case "character":
                return f"character(len=1, kind=c_char)"
            case "integer":
                return f"integer(c_int)"
            case "real":
                return f"real(c_double)"
            case "complex":
                return f"real(c_double_complex)"
            case _:
                raise ValueError(f"'{self.name}' not supported for C formatting")

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
            return f"dimension({", ".join(self)}), "

    def Intent(self, spec):
        return f"intent({self.name.lower()})"

    def ProcedureArgument(self, spec):
        pass

    def Procedure(self, spec):
        match spec:
            case "call":
                call = f"{self.name}({", ".join(arg.name for arg in self.args)})"
                if self.retvar is None:
                    call = "call " + call
                else:
                    call = f"{self.retvar} = {call}"
                return Indentable(call)
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

    def C_Argument(self, spec):
        match spec:
            case "dummy":
                arg_str = f"{format(self.type, "C")}, {self.dimension}{self.intent}, target :: {self.name}"
                for line in self.doc:
                    arg_str += f"\n    !!{line}"
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")
        return Indentable(arg_str)

    def C_Arguments(self, spec):
        match spec:
            case "dummy":
                formatted = "\n".join(format(arg, "dummy") for arg in self)
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
        doc = [f" C-wrapper for {self.orig_proc_ford_link}"] + self.doc

        wrapper = f"""{format(doc, "subroutine")}
subroutine {self.name}({format(self.arguments, "arglist")}) bind(C, name="{self.name}")
    use {self.module_name}, only: {self.orig_proc_name}
{format(self.arguments, "dummy") >> INDENT}

{format(self.arguments, "null_validation") >> INDENT}

{self.orig_proc_call >> INDENT}
end subroutine {self.name}"""

        return Indentable(wrapper)

    def C_Module(self, spec):
        doc = [f" Module for C-wrappers for [[{self.orig_mod_name}(module)]]"] + self.doc
        return f"""#ifndef NO_C_INTERFACE
#include <src/macros.h>

{format(doc, "module")}
module {self.name}
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

{"\n\n    ".join(format(c_wrapper) >> INDENT for c_wrapper in self)}

end module {self.name}
#endif"""
