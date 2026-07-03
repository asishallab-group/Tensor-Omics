import warnings


def eval_expr(expr: str):
    expr = expr.lower()

    expr = expr.replace("acos", "np.arccos")
    expr = expr.replace("acosh", "np.arccosh")
    expr = expr.replace("asin", "np.arcsin")
    expr = expr.replace("asinh", "np.arcsinh")
    expr = expr.replace("atan", "np.arctan")
    expr = expr.replace("atan2", "np.arctan2")
    expr = expr.replace("atanh", "np.arctanh")
    expr = expr.replace("cos", "np.cos")
    expr = expr.replace("cosh", "np.cosh")
    expr = expr.replace("sin", "np.sin")
    expr = expr.replace("sinh", "np.sinh")
    expr = expr.replace("tan", "np.tan")
    expr = expr.replace("tanh", "np.tanh")
    expr = expr.replace("achar", "chr")
    expr = expr.replace("char", "chr")

    expr = expr.replace("_int32", "")
    expr = expr.replace("_real64", "")

    expr = expr.replace(".true.", "")
    expr = expr.replace(".false.", "")

    expr = expr.replace("pi", str(np.pi))

    return eval(expr)


def warn(msg, entity):
    warnings.warn(extend_err_msg(f"\n\033[38;5;226m{msg}", entity))


def error(cls, msg, entity):
    raise cls(extend_err_msg(f"\n\033[38;5;196m{msg}", entity))


def extend_err_msg(msg, entity):
    extended = f"{msg}\033[38;5;226m for '\033[38;5;210m{entity.name}\033[38;5;226m'"
    for parent in iter_parents(entity):
        extended += f" in '\033[38;5;208m{parent.name}\033[38;5;226m'"
    return extended + "\033[0m"


def iter_parents(entity):
    parent = entity
    while (parent := getattr(parent, "parent", None)) is not None:
        yield parent


def repeat(value, n_times=None):
    if type(n_times) is int:
        for i in range(n_times):
            yield value
    else:
        while True:
            yield value


class Indentable(str):
    """Wrapper class to easily indent code blocks"""
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *args, **kwargs)

    def __rshift__(self, level):
        indented = "\n".join(level * " " + line  if line else "" for line in self.split("\n"))
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

    def FordLink(self, spec):
        return str(self)

    def TextAlignment(self, length):
        length = int(length)
        match self:
            case self.CENTER:
                return ":" + "-" * (length - 2) + ":"
            case self.RIGHT:
                return "-" * (length - 1) + ":"
            case self.LEFT:
                return "-" * length

    def FordTable(self, spec):
        header = tuple(format(col, spec) for col in self.header)
        rows = tuple(tuple(format(col, spec) for col in row) for row in self.rows)
        col_widths = tuple(max(len(header_col), max(map(len, rows_col))) for header_col, rows_col in zip(header, zip(*rows)))

        header = " | ".join(col.center(width, " ") for col, width in zip(header, col_widths))
        align = "|".join(format(a, str(width + 2)) for a, width in zip(self.align, col_widths))
        rows = tuple(" | ".join(col.center(width, " ") for col, width in zip(row, col_widths)) for row in rows)

        header = f"| {header} |"
        align = f"|{align}|"

        return Indentable("\n".join((header, align, *(f"| {row} |" for row in rows))))

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
            serializer_func = getattr(serializer, cls.__name__)
            try:
                return serializer_func(self, spec)
            except NotImplementedError as e:
                if not str(e):
                    raise NotImplementedError(f"No implementation for Serializer.{cls.__name__}(self, spec)")
                else:
                    raise e


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
