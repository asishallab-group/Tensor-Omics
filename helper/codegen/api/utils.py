class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super(cls, cls).__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        return level * " " + self.replace("\n", "\n" + level * " ")


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
