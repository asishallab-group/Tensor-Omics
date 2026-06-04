from ford.fortran_project import Project as _Project
from ford.settings import load_toml_settings, load_markdown_settings
from ford.sourceform import FortranSubroutine, FortranFunction, FortranVariable, FortranModuleProcedureImplementation
import warnings
from enum import Enum
import re
from .utils import Indentable


class Project(_Project):
    """Class for the parsed codebase"""
    def __init__(self):
        directory = "."

        # load settings from fpm.toml
        proj_settings = load_toml_settings(directory)

        # if no fpm.toml, use ford.yml
        if proj_settings is None:
            with open("ford.yml", "r") as f:
                proj_settings, _ = load_markdown_settings(directory, f.read(), f.name)

        super(Project, self).__init__(proj_settings)
        self.correlate()


# TODO: handle derived types
class Fortran_Type:
    def __init__(self, type_name: str, kind: str = None, length: str = None):
        if type_name not in ("character", "logical", "type") and kind is None:
            raise ValueError(f"kind must be present for '{type_name}'")
        self.name = type_name
        self.kind = kind
        self.len = length

    @classmethod
    def from_fortran_variable(cls, arg: FortranVariable):
        type = re.match(r"\s*([^ \(]+).*", arg.full_type, re.IGNORECASE).group(1).lower()

        if type == "character":
            len_match = re.match(r".*\blen\s*=\s*(\*|[^\),]+).*", arg.full_type, re.IGNORECASE)
            if len_match is None:
                raise SyntaxError("Dummy arguments of type 'character' need to be declared with a named 'len' type paramter, like 'character(len=42)'")
            length = len_match.group(1)
        else:
            length = None

        try:
            return cls(type, arg.kind, length)
        except ValueError as e:
            raise ValueError(f"{e} of '{arg.name}' in '{arg.parent.name}'")

    def __format__(self, spec):
        match spec:
            case "C":
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
            case "Fortran":
                match self.name:
                    case "character":
                        return f"character(len={self.length})"
                    case "logical":
                        return "logical"
                    case "type":
                        raise ValueError(f"'{self.name}' not supported for C formatting")
                    case _:
                        return f"{self.type}(kind={self.kind})"
            case _:
                raise ValueError(f"Unsupported formatting spec '{spec}'")


class DocList:
    """Class for managing parsed Ford documentation. Each line is one element"""
    def __init__(self, doc_list: List[str, ...], type: str):
        if type not in ("module", "procedure", "argument"):
            raise ValueError("type must be one of: 'module', 'procedure', 'argument'")
        self.doc_list = doc_list
        self.type = type

    @classmethod
    def from_fortran(cls, unit: FortranModule | FortranSubroutine | FortranFunction | FortranVariable | FortranModuleProcedureImplementation):
        match type(unit).__name__:
            case "FortranModule":
                ty = "module"
            case "FortranSubroutine" | "FortranFunction":
                ty = "procedure"
            case "FortranVariable":
                ty = "argument"
            case _:
                raise TypeError(f"DocList doesn't support '{type(unit).__name__}'")

        doc_list = unit.doc_list
        if len(doc_list) == 0 or (len(doc_list) == 1 and doc_list[0] == ""):
            warnings.warn(f"\n\033[38;5;226mNo docstring for '\033[38;5;208m{unit.name}\033[38;5;226m' in '{unit.parent.name}\033[0m'")

        return cls(doc_list, ty)

    def __getitem__(self, idx):
        return self.doc_list[idx]

    def __setitem__(self, idx, value):
        self.doc_list[idx] = value

    def __str__(self):
        match self.type:
            case "argument":
                formatted = f"!!{"\n!!".join(self.doc_list)}"
            case "procedure" | "module":
                formatted = f"!|{"\n!|".join(self.doc_list)}"
        return Indentable(arg_str)

    def __add__(self, other: list):
        return DocList([*new.doc_list, *other], self.type)

    def __radd__(self, other: list):
        return DocList([*other, *new.doc_list], self.type)


class Dimension(tuple):
    """Representation of the dimension attribute of a variable"""
    def __new__(cls, arg: List[str, ...] = ()):
        return super(cls, cls).__new__(cls, arg)

    @classmethod
    def from_fortran_variable(cls, arg: FortranVariable):
        dims = arg.dimension
        for attrib in arg.attribs:
            if (match := re.match(r".*dimension\s*(?P<dimension>\([^\)]+\)).*", attrib)) is not None:
                dims = match.group("dimension")

        dims = dims[1:-1]

        if not dims:
            return cls()
        else:
            # TODO: shape might include some funtion like 'int(n, int32)' -> results in wrong extent
            return cls([dim.strip() for dim in dims.split(",")])

    def __str__(self):
        if len(self) == 0:
            return ""
        else:
            return f"dimension({", ".join(self)}), "


class Intent(Enum):
    IN = 1
    OUT = 2
    INOUT = 3

    @classmethod
    def _missing_(cls, value):
        if type(value) is FortranVariable:
            intent = value.intent.upper()
            for member in cls:
                if member.name == intent:
                    return member
            raise SyntaxError(f"No intent for '{value.name}' in '{value.parent.name}' of module '{value.parent.parent.name}'")
        raise KeyError(f"'{value}' invalid intent")

    def __str__(self):
        return f"intent({self.name.lower()})"


class ProcedureArgument:
    """Wrapper class for procedure arguments"""

    def __init__(self, argument: FortranVariable):
        self.name = argument.name
        self.doc_list = DocList.from_fortran(argument)
        self.intent = Intent(argument)
        self.dimension = Dimension.from_fortran_variable(argument)
        self.attribs = argument.attribs
        self.full_type = argument.full_type
        self.optional = argument.optional
        self.type = Fortran_Type.from_fortran_variable(argument)

    @property
    def default_value(self) -> int | str | bool | float | None:
        # TODO: parse default value from doc_list[-1]
        if self.optional:
            ...


class Procedure:
    """Wrapper class for Fortran procedures"""
    def __init__(self, procedure: FortranSubroutine | FortranFunction | FortranModuleProcedureImplementation):
        self.name = procedure.name
        self.meta = procedure.meta
        self.args = tuple(map(ProcedureArgument, procedure.args))
        self.doc_list = DocList.from_fortran(procedure)
        self.retvar = getattr(procedure, "retvar", None)
        self.type = "subroutine" if self.retvar is None else "function"
        self.parent = procedure.parent

    def __format__(self, spec):
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
