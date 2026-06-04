from abc import ABC, abstractmethod


class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super(cls, cls).__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        return level * " " + self.replace("\n", "\n" + level * " ")


class Serializer(ABC):
    @classmethod
    @abstractmethod
    def Fortran_Type(self, spec):
        pass

    @classmethod
    @abstractmethod
    def DocList(self, spec):
        pass

    @classmethod
    @abstractmethod
    def Dimension(self, spec):
        pass

    @classmethod
    @abstractmethod
    def Intent(self, spec):
        pass

    @classmethod
    @abstractmethod
    def ProcedureArgument(self, spec):
        pass

    @classmethod
    @abstractmethod
    def Procedure(self, spec):
        pass

    @classmethod
    @abstractmethod
    def C_Argument(self, spec):
        pass

    @classmethod
    @abstractmethod
    def C_Arguments(self, spec):
        pass

    @classmethod
    @abstractmethod
    def C_Wrapper(self, spec):
        pass

    @classmethod
    @abstractmethod
    def C_Module(self, spec):
        pass

    @classmethod
    @abstractmethod
    def C_Modules(self, spec):
        pass


class CodeGenerator:
    @classmethod
    def use(cls, serializer: Serializer):
        if serializer.__bases__[0] is not Serializer:
            raise TypeError(f"'serializer' must have 'Serializer' as first base class")
        CodeGenerator.serializer = serializer

    def __format__(self, spec):
        cls = type(self)
        serialier_func = getattr(CodeGenerator.serializer, cls.__name__)
        return serialier_func(self, spec)
