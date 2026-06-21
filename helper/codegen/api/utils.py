import warnings


def warn(msg, entity):
    warning = f"\n\033[38;5;226m{msg} for '\033[38;5;210m{entity.name}\033[38;5;226m'"
    parent = entity
    while (parent := getattr(parent, "parent", None)) is not None:
        warning += f" in '\033[38;5;208m{parent.name}\033[38;5;226m'"

    warnings.warn(warning + "\033[0m")


class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        indented = "\n".join(level * " " + line for line in self.split("\n") if line)
        return Indentable(indented)


class Snippet:
    """Class for VS Code compatible Snippet"""
    def __init__(self, key, prefix, body, description):
        self.key = key
        self.prefix = prefix
        self.body = body
        self.description = description

    @classmethod
    def create_setup_snippet(cls, language, language_short, snippet_type, routine_name, module_name, body, description):
        return cls(
            key=f"{language} setup for '{routine_name}' from module '{module_name}'",
            prefix=f"{snippet_type}:{language_short}_setup_{routine_name}",
            body=body.split("\n"),
            description=description
        )

    @classmethod
    def create_call_snippet(cls, language, language_short, snippet_type, routine_name, module_name, body, description):
        return cls(
            key=f"{language} call of '{routine_name}' from module '{module_name}'",
            prefix=f"{snippet_type}:{language_short}_{routine_name}",
            body=body.split("\n"),
            description=description
        )


class Serializer:
    @classmethod
    def dump(self, *args, **kwargs):
        raise NotImplementedError()

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

    def Error_Handling(self, spec):
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
        self.preprocessor = CmdPreprocessor(["", macro_file])

        if regex_escaped:
            from re import escape
            for macro in self.preprocessor.macros.values():
                for token in macro.value:
                    token.value = escape(token.value)

    def expand(self, text: str):
        tokens = self.preprocessor.tokenize(text)
        expanded = self.preprocessor.expand_macros(tokens)
        return "".join(i.value for i in expanded)


regex_escaped_preprocessor = Macros("src/macros.h")
preprocessor = Macros("src/macros.h", False)
