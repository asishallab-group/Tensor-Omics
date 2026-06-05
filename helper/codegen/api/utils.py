class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super(cls, cls).__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        indented = ""
        for line in self.split("\n"):
            if line:
                indented += level * " " + line
            indented += "\n"
        return Indentable(indented)


class Serializer:
    def Module(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Module(self, spec):
        raise NotImplementedError()

    def Modules(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Modules(self, spec):
        raise NotImplementedError()

    def Fortran_Type(self, spec):
        raise NotImplementedError()

    def DocList(self, spec):
        raise NotImplementedError()

    def Dimension(self, spec):
        raise NotImplementedError()

    def Intent(self, spec):
        raise NotImplementedError()

    def Procedure(self, spec):
        raise NotImplementedError()

    def C_Wrapper(self, spec):
        raise NotImplementedError()

    def Procedure_Argument(self, spec):
        raise NotImplementedError()

    def Procedure_Arguments(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Argument(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Arguments(self, spec):
        raise NotImplementedError()


class CodeGenerator:
    @classmethod
    def use(cls, serializer: Serializer):
        if serializer.__bases__[0] is not Serializer:
            raise TypeError(f"'serializer' must have 'Serializer' as first base class")
        CodeGenerator.serializer = serializer

    def __format__(self, spec):
        if (serializer := getattr(CodeGenerator, "serializer", None)) is None:
            return Indentable("")
        else:
            cls = type(self)
            serialier_func = getattr(serializer, cls.__name__)
            try:
                return serialier_func(self, spec)
            except NotImplementedError as e:
                raise NotImplementedError(f"No implementation for Serializer.{cls.__name__}(self, spec)")


class Macros:
    """Class for Macros"""
    def __init__(self, macro_file: str, regex_escaped=True):
        from pcpp.pcmd import CmdPreprocessor
        from re import escape

        self.preprocessor = CmdPreprocessor(["", macro_file])
        if regex_escaped:
            for macro in self.preprocessor.macros.values():
                for token in macro.value:
                    token.value = escape(token.value)

    def __format__(self, spec):
        tokens = self.preprocessor.tokenize(spec)
        expanded = self.preprocessor.expand_macros(tokens)
        return "".join(i.value for i in expanded)


regex_escaped_preprocessor = Macros("src/macros.h")
preprocessor = Macros("src/macros.h", False)
