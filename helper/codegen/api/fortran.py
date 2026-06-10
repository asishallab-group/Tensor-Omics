from ford.fortran_project import Project as _Project
from ford.settings import load_toml_settings, load_markdown_settings
from ford.sourceform import FortranSubroutine, FortranFunction, FortranVariable, FortranModuleProcedureImplementation, FortranModule
import warnings
from enum import Enum
import re
from .utils import Indentable, CodeGenerator, regex_escaped_preprocessor


SHAPE_ARG_SUFFIX = "_shape"


class Module(CodeGenerator):
    """Class for a parsed Fortran Module"""
    def __init__(self, module: FortranModule):
        self.name = module.name
        self.doc_list = DocList.from_fortran(module)
        self.procedures = tuple(map(Procedure, module.routines))


class Modules(CodeGenerator, tuple):
    """Collection of Module objects"""
    def __new__(cls):
        proj = Project()

        # sort modules for deterministic order
        sorted_mods = sorted(proj.modules, key=lambda x: x.name)

        mods = super().__new__(cls, map(Module, sorted_mods))
        mods.project = proj
        return mods


class Project(CodeGenerator, _Project):
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


class Dimension(CodeGenerator, tuple):
    """Representation of the dimension attribute of a variable"""
    def __new__(cls, arg: List[str, ...] = ()):
        return super().__new__(cls, arg)

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

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if isinstance(item, slice):
            return Dimension(result)
        return result


class Intent(CodeGenerator, Enum):
    IN = 1
    OUT = 2
    INOUT = 3

    @classmethod
    def _missing_(cls, value):
        if type(value) is FortranVariable:
            intent = value.intent.upper()
            if intent:
                for member in cls:
                    if member.name == intent:
                        return member
            else:
                return cls.INOUT
            raise SyntaxError(f"No intent for '{value.name}' in '{value.parent.name}' of module '{value.parent.parent.name}'")
        raise KeyError(f"'{value}' invalid intent")


# TODO: handle derived types
class Fortran_Type(CodeGenerator):
    def __init__(self, type_name: str, intent: Intent, kind: str = None, dimension: Dimension = Dimension(), length: str = None):
        if type_name not in ("character", "logical", "type") and kind is None:
            raise ValueError(f"kind must be present for '{type_name}'")
        self.name = type_name
        self.kind = kind
        self.len = length
        self.is_fixed_length = False if length is None else re.match(r"\d+", length)
        self.dimension = dimension
        self.intent = intent
        self.needs_conversion = self.name in ("character", "logical")

    @classmethod
    def from_fortran_variable(cls, arg: FortranVariable):
        type = re.match(r"\s*([^ \(]+).*", arg.full_type, re.IGNORECASE).group(1).lower()

        if type == "character":
            len_match = re.match(r".*\blen\s*=\s*(\*|[^\), ]+).*", arg.full_type, re.IGNORECASE)
            if len_match is None:
                raise SyntaxError("Dummy arguments of type 'character' need to be declared with a named 'len' type paramter, like 'character(len=42)'")
            length = len_match.group(1)
        else:
            length = None

        try:
            return cls(type, intent=Intent(arg), kind=arg.kind, dimension=Dimension.from_fortran_variable(arg), length=length)
        except ValueError as e:
            raise ValueError(f"{e} of '{arg.name}' in '{arg.parent.name}' in '{arg.parent.parent.name}'")


class DocList(CodeGenerator):
    """Class for managing parsed Ford documentation. Each line is one element"""
    def __init__(self, doc_list: List[str, ...], type: str):
        if type not in ("module", "procedure", "argument"):
            raise ValueError("type must be one of: 'module', 'procedure', 'argument'")
        self.doc_list = list(doc_list)
        self.type = type

        while len(self.doc_list) > 0 and self.doc_list[-1] == "":
            self.doc_list.pop()

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
            warning = f"\n\033[38;5;226mNo docstring for '\033[38;5;208m{unit.name}\033[38;5;226m'"
            parent = unit
            while type(parent) is not FortranModule:
                parent = parent.parent
                warning += f" in '\033[38;5;208m{parent.name}\033[38;5;226m'"

            warnings.warn(warning + "\033[0m")

        doc_list = [line.lstrip() for line in doc_list]
        return cls(doc_list, ty)

    def __getitem__(self, idx):
        return self.doc_list[idx]

    def __setitem__(self, idx, value):
        self.doc_list[idx] = value

    def __add__(self, other: list):
        return DocList([*self.doc_list, *other], self.type)

    def __radd__(self, other: list):
        return DocList([*other, *self.doc_list], self.type)

    def __len__(self):
        return len(self.doc_list)


class Procedure_Argument(CodeGenerator):
    """Wrapper class for procedure arguments"""

    def __init__(self, argument: FortranVariable, proc: Procedure):
        self.name = argument.name
        self.doc_list = DocList.from_fortran(argument)
        self.attribs = argument.attribs
        self.optional = argument.optional
        self.type = Fortran_Type.from_fortran_variable(argument)
        self.is_temporary = self.name.startswith("tmp_")
        self.parent = proc
        self.is_shape_arg = self.name.endswith(SHAPE_ARG_SUFFIX)

    @property
    def default_value(self) -> int | str | bool | float | None:
        if len(self.doc_list) > 0:
            doc_str = self.doc_list[-1]
            regex = format(regex_escaped_preprocessor, ".*M_DOC_DEFAULT((?P<default_val>.*)).*")
            if (match := re.match(regex, doc_str)) is not None:
                return match.group("default_val")

        return None

    @property
    def is_dim_arg_for(self) -> Self | None:
        if len(self.type.dimension) == 0:
            for arg in self.parent.args:
                if self.name in arg.type.dimension:
                    return arg

    @property
    def shape_arg(self) -> Self | None:
        if self.type.dimension == Dimension((":")):
            shape_arg_name = self.name + SHAPE_ARG_SUFFIX
            for arg in self.parent.args:
                if arg.name == shape_arg_name:
                    if len(arg.type.dimension) == 0:
                        raise SyntaxError(f"Shape argument '{self.name}' must have a dimension attribute in {arg.parent.name} of {arg.parent.parent.name}")
                    if arg.type.intent is not Intent.IN:
                        raise SyntaxError(f"Shape argument '{self.name}' must have 'intent(in)' in {arg.parent.name} of {arg.parent.parent.name}")
                    return arg


class Procedure_Arguments(CodeGenerator, tuple):
    """Collection of Procedure_Argument objects"""
    def __new__(cls, args: List[FortranVariable], proc: Procedure):
        return super().__new__(cls, (Procedure_Argument(arg, proc) for arg in args))


class Procedure(CodeGenerator):
    """Wrapper class for Fortran procedures"""
    def __init__(self, procedure: FortranSubroutine | FortranFunction | FortranModuleProcedureImplementation):
        self.name = procedure.name
        self.meta = procedure.meta
        self.args = Procedure_Arguments(procedure.args, self)
        self.doc_list = DocList.from_fortran(procedure)
        self.retvar = getattr(procedure, "retvar", None)
        self.type = "subroutine" if self.retvar is None else "function"
        self.parent = procedure.parent
