from ford.fortran_project import Project
from ford.sourceform import FortranSubroutine, FortranFunction, FortranVariable
from typing import List, Tuple
from enum import Enum


class C_Type(Enum):
    C_INT = 1
    C_DOUBLE = 2
    C_DOUBLE_COMPLEX = 3
    C_CHAR = 4

    @classmethod
    def _missing_(cls, value):
        match value:
            case "real64": return cls.C_DOUBLE
            case "int32": return cls.C_INT
            case "complex(real64)": return cls.C_DOUBLE_COMPLEX
            case "character": return cls.C_CHAR
            case _:
                raise KeyError(f"'{value}' not valid for c types")


class Intent(Enum):
    IN = 1
    OUT = 2
    INOUT = 3

    @classmethod
    def _missing_(cls, value):
        value = value.upper()
        for member in cls:
            if member.value == value:
                return member
        raise KeyError(f"'{value}' invalid intent")


Dimension = Tuple[str, ...]


class C_Argument:
    """Includes all relevant information about an argument, meaning its type, dimension, name, docstring and default value if optional"""
    def __init__(self, name: str, doc: List[str], type: C_Type, dimension: Dimension, intent: Intent, default_value=None):
        self.name = name
        self.doc = doc
        self.type = type
        self.dimension = dimension
        self.intent = intent
        self.default_value = default_value


C_Arguments = Tuple[C_Argument]


class C_Wrapper:
    """Includes all relevant information about a C wrapper, meaning its C name, arguments and argument docstrings"""
    def __init__(self, name: str, doc: List[str], arguments: C_Arguments):
        self.name = name
        self.doc = doc
        self.arguments = arguments


def generate_c_wrappers(proj: Project, out_dir=None) -> List[C_Wrapper]:
    """
        Generates the C wrappers for the subroutines/functions in the Project that have the 'category: C-interface' meta-tag
        and writes them to out_dir if present (else only returns Wrappers)
    """

    name_suffix = "_c"

    for module in proj.modules:
        for subroutine in module.subroutines:
            if subroutine.meta.category == "C-interface":
                name = get_wrapper_name_subroutine(subroutine, name_suffix)
                doc = subroutine.doc_list
                arguments = get_wrapper_arguments_subroutine(subroutine)


def get_wrapper_arguments_subroutine(subroutine: FortranSubroutine) -> C_Arguments:
    return get_wrapper_arguments(subroutine.args)


def get_wrapper_arguments_function(function: FortranFunction) -> C_Arguments:
    arg = function.retvar
    result_argument = C_Argument(
        name=arg.name,
        doc=arg.doc_list,
        type=get_c_type(arg),
        dimension=get_dimension(arg),
        intent=Intent.OUT
    )

    return tuple(...get_wrapper_arguments(function.args), result_argument)


def get_wrapper_arguments(args: List[FortranVariable]) -> C_Arguments:
    arguments = []
    for arg in args:
        arguments.append(C_Argument(
            name=arg.name,
            doc=arg.doc_list,
            type=get_c_type(arg),
            dimension=get_dimension(arg),
            intent=Intent(arg.intent),
            default_value=get_default_value(arg)
        ))

    return tuple(arguments)


def get_default_value(arg: FortranVariable) -> int | str | bool | float | None:
    # TODO: parse default value from doc_list[-1]
    if arg.optional:
        ...


def get_dimension(arg: FortranVariable) -> Dimension:
    dims = arg.dimension[1:-1]
    dimension = []
    for dim_idx, dim in enumerate(dims.split(","), 1):
        dim = dim.strip()
        if dim == ":":
            dimension.append(f"n_{arg.name}_elements_dim_{dim_idx}")
        else:
            dimension.append(dim)

    return tuple(dimension)


def get_c_type(arg: FortranVariable) -> str:
    if "complex(" in arg.full_type:
        type_name = f"complex({arg.kind})"
    elif "character" in arg.full_type:
        type_name = "character"
    else:
        type_name = arg.kind
    return C_Type(type_name)


def get_wrapper_name_subroutine(subroutine: FortranSubroutine, name_suffix: str) -> str:
    """
        Creates the name for the C wrapper of the passed subroutine, respecting the expert/non-expert convention
    """
    # alloc routines don't get the alloc suffix in the C wrapper name
    if subroutine.name.endswith("_alloc"):
        name = subroutine.name.removesuffix("_alloc")
    # non-alloc routines get the _expert suffix if there is an alloc variant
    elif module.find_child(f"{subroutine.name}_alloc") is not None:
        name = f"{subroutine.name}_expert"
    # if no alloc variant, the base name is simply the subroutine name
    else:
        name = subroutine.name

    return name + name_suffix
