from ford.fortran_project import Project as _Project
from ford.settings import load_toml_settings, load_markdown_settings
from ford.sourceform import FortranSubroutine, FortranFunction, FortranVariable, FortranModuleProcedureImplementation, FortranModule
from enum import Enum
import re
from .utils import Indentable, CodeGenerator, warn, error, repeat
from .doc import DocList, FordTable
import numpy as np


SHAPE_ARG_SUFFIX = "_shape"


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


class Error_Handling(CodeGenerator):
    """Class for everything related to error handling"""
    def __init__(self, project: Project):
        tox_errors = project.find("tox_errors")
        self.error_codes = {
            var.name: (int(eval_expr(var.initial)), DocList(var, var.doc_list, "variable"))
            for var in tox_errors.variables if var.parameter and var.name.startswith("ERR_")
        }


class Module(CodeGenerator):
    """Class for a parsed Fortran Module"""
    def __init__(self, module: FortranModule, error_handling: Error_Handling):
        self.name = module.name
        self.doc_list = DocList.from_fortran(self, module)
        self.procedures = tuple(map(Procedure, module.routines, repeat(self)))
        self.error_handling = error_handling
        self.parent = module.parent


class Modules(CodeGenerator, tuple):
    """Collection of Module objects"""
    def __new__(cls, exclude_directories=[]):
        proj = Project(exclude_directories)

        # sort modules for deterministic order
        sorted_mods = sorted(proj.modules, key=lambda x: x.name)

        error_handling = Error_Handling(proj)

        mods = super(Modules, cls).__new__(cls, map(Module, sorted_mods, repeat(error_handling)))
        mods.project = proj

        mods.error_handling = error_handling
        return mods


class Project(CodeGenerator, _Project):
    """Class for the parsed codebase"""
    def __init__(self, exclude_directories=[]):
        directory = "."

        # load settings from fpm.toml
        proj_settings = load_toml_settings(directory)

        # if no fpm.toml, use ford.yml
        if proj_settings is None:
            with open("ford.yml", "r") as f:
                proj_settings, _ = load_markdown_settings(directory, f.read(), f.name)

        proj_settings.exclude_dir.extend(exclude_directories)

        super(Project, self).__init__(proj_settings)

        self.correlate()


class Dimension(CodeGenerator, tuple):
    """Representation of the dimension attribute of a variable"""
    def __new__(cls, arg: List[str, ...] = ()):
        return super(Dimension, cls).__new__(cls, arg)

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

    @property
    def ndim(self):
        return len(self.dimension)

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


class Procedure_Argument(CodeGenerator):
    """Wrapper class for procedure arguments"""

    # e.g. mode or mode_<arg>, or method or method_<arg>
    MODE_ARG_RE = re.compile(rf"^(?P<alias>{"|".join(map(str.lower, FordTable.MODE_ALIASES))})(?:_.+)?$")

    def __init__(self, argument: FortranVariable, proc: Procedure):
        self.parent = proc
        self.name = argument.name
        self.doc_list = DocList.from_fortran(self, argument)
        self.attribs = argument.attribs
        self.optional = argument.optional
        self.type = Fortran_Type.from_fortran_variable(argument)
        self.is_temporary = self.name.startswith("tmp_")
        self.is_shape_arg = self.name.endswith(SHAPE_ARG_SUFFIX)

        mode_arg_match = self.MODE_ARG_RE.match(self.name)
        self.is_mode_arg = mode_arg_match is not None
        if self.is_mode_arg:
            if self.doc_list.meta["mode_table"] is None:
                alias = mode_arg_match.group("alias").capitalize()
                error(SyntaxError, FordTable.mode_var_table_example(alias), self)

        # determine default value
        self.default_value = self.doc_list.meta["default"]
        if self.default_value is not None:
            self.default_value = eval_expr(self.default_value.group("default_val"))

    @property
    def is_mask_count_arg_for(self) -> Self | None:
        match = re.match(r"n_selected_(?P<mask_arg_name>.*)", self.name)
        if match is not None:
            mask_arg_name = match.group("mask_arg_name")
            for arg in self.parent.args:
                if re.match(f"{mask_arg_name}_(?:selection_)?mask", arg.name) is not None:
                    return arg

    @property
    def is_dim_arg_for(self) -> Tuple[Self, ...]:
        dim_args = []
        if len(self.type.dimension) == 0:
            for arg in self.parent.args:
                if self.name in arg.type.dimension:
                    dim_args.append(arg)
        return tuple(dim_args)

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
        return super(Procedure_Arguments, cls).__new__(cls, map(Procedure_Argument, args, repeat(proc)))


class Procedure(CodeGenerator):
    """Wrapper class for Fortran procedures"""
    def __init__(self, procedure: FortranSubroutine | FortranFunction | FortranModuleProcedureImplementation, parent: Module):
        self.parent = parent
        self.fortran_procedure = procedure
        self.name = procedure.name
        self.meta = procedure.meta
        if not self.meta.summary:
            warn("No summary meta tag (!! summary: ...)", procedure)
        if not self.meta.author:
            warn("No author meta tag (!! author: ...)", procedure)
        self.args = Procedure_Arguments(procedure.args, self)
        self.doc_list = DocList.from_fortran(self, procedure)
        self.retvar = getattr(procedure, "retvar", None)
        if type(self.retvar) is FortranVariable:
            self.retvar = Procedure_Argument(self.retvar, self)
            self.retvar.type.intent = Intent.OUT
        self.type = "subroutine" if self.retvar is None else "function"
