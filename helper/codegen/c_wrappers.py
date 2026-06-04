from ford.fortran_project import Project
from ford.sourceform import FortranSubroutine, FortranFunction, FortranVariable, FortranModuleProcedureImplementation
from typing import List, Tuple
from enum import Enum
import re
from pathlib import Path
import os


Procedure = FortranSubroutine | FortranFunction | FortranModuleProcedureImplementation


class C_Type(Enum):
    C_INT = 1
    C_DOUBLE = 2
    C_DOUBLE_COMPLEX = 3
    C_CHAR = 4
    C_BOOL = 5

    @classmethod
    def _missing_(cls, value):
        match value:
            case "real64": return cls.C_DOUBLE
            case "int32": return cls.C_INT
            case "complex(real64)": return cls.C_DOUBLE_COMPLEX
            case "character": return cls.C_CHAR
            case _:
                raise KeyError(f"'{value}' not valid for c types")

    def __str__(self):
        match self:
            case self.C_BOOL:
                return f"integer(c_int)"
            case self.C_CHAR:
                return f"character(len=1, kind=c_char)"
            case self.C_INT:
                return f"integer(c_int)"
            case self.C_DOUBLE:
                return f"real(c_double)"
            case self.C_DOUBLE_COMPLEX:
                return f"real(c_double_complex)"


class Intent(Enum):
    IN = 1
    OUT = 2
    INOUT = 3

    @classmethod
    def _missing_(cls, value):
        value = value.upper()
        for member in cls:
            if member.name == value:
                return member
        raise KeyError(f"'{value}' invalid intent")

    def __str__(self):
        return f"intent({self.name.lower()})"


Dimension = Tuple[str, ...]


class C_Argument:
    """Includes all relevant information about an argument, meaning its type, dimension, name, docstring and default value if optional"""
    def __init__(self, name: str, doc: List[str], type: C_Type, dimension: Dimension, intent: Intent, default_value=None, only_c_arg=False):
        self.name = name
        self.doc = doc
        self.type = type
        self.dimension = dimension
        self.intent = intent
        self.default_value = default_value
        self.only_c_arg = only_c_arg

    def __str__(self):
        dimension = "" if len(self.dimension) == 0 else f"dimension({", ".join(self.dimension)}), "
        return f"""{self.type}, {dimension}{self.intent}, target :: {self.name}
    !!{"\n    !!".join(self.doc)}"""


C_Arguments = Tuple[C_Argument]


class C_Wrapper:
    """Includes all relevant information about a C wrapper, meaning its C name, arguments and argument docstrings"""
    def __init__(self, name: str, orig_name: str, doc: List[str], arguments: C_Arguments, module_name: str, retvar: str | None):
        self.name = name
        self.orig_name = orig_name
        self.doc = doc
        self.arguments = arguments
        self.module_name = module_name
        self.retvar = retvar

    def __str__(self):

        call = f"{self.orig_name}({", ".join(arg.name for arg in self.arguments if not arg.only_c_arg)})"
        if self.retvar is None:
            call = "call " + call
        else:
            call = f"{self.retvar} = {call}"

        return f"""!>{"\n!|".join(self.doc)}
subroutine {self.name}({", ".join(arg.name for arg in self.arguments)}) bind(C, name="{self.name}")
    use {self.module_name}, only: {self.orig_name}
{indent("\n".join(str(arg) for arg in self.arguments), 4)}

    M_CHECK_IERR_NON_NULL
    {"\n    ".join(f"M_CHECK_NON_NULL({arg.name})" for arg in self.arguments if arg.name != "ierr")}

    {call}
end subroutine {self.name}"""


class C_Module:
    """Includes a bunch of wrappers that should appear in one c interfacing module"""
    def __init__(self, name: str, doc: List[str]):
        self.name = name
        self.doc = doc
        self.c_wrappers = []

    def __iadd__(self, c_wrapper):
        if type(c_wrapper) is C_Wrapper:
            self.c_wrappers.append(c_wrapper)
        else:
            raise TypeError("Only supporting values of type 'C_Wrapper' to be added")

        return self

    def __iter__(self):
        return iter(self.c_wrappers)

    def __len__(self):
        return len(self.c_wrappers)

    def __str__(self):
        return f"""#include <src/macros.h>

!>{"\n!|".join(self.doc)}
module {self.name}
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

{"\n\n    ".join(indent(str(c_wrapper), 4) for c_wrapper in self)}

end module {self.name}"""


C_Modules = List[C_Module]


def indent(code: str, level=0) -> str:
    return level * " " + code.replace("\n", "\n" + level * " ")


def generate_c_module_code(c_modules: C_Modules, out_dir: str):
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise ValueError("out_dir muste be a valid directory path")

    c_wrapper_dir = out_dir.joinpath("c_interface")
    if c_wrapper_dir.is_dir():
        from shutil import rmtree
        rmtree(c_wrapper_dir)

    os.mkdir(c_wrapper_dir)

    for module in c_modules:
        module_file_path = c_wrapper_dir.joinpath(module.name + ".F90")

        with open(module_file_path, "w") as module_file:
            module_file.write(str(module))


def collect_c_modules(project: Project) -> C_Modules:
    """
        Creates and Collects all C wrappers for the subroutines/functions in the Project that have the 'category: C-interface' meta-tag and returns them.
        The output will only have non-empty modules.
    """

    name_suffix = "_c"

    c_modules = []

    for module in project.modules:
        if not module.name.endswith(name_suffix):
            c_module = C_Module(module.name + name_suffix, module.doc_list)
            for procedure in module.routines:
                if procedure.meta.category == "C-interface":
                    name = get_wrapper_name(procedure, name_suffix)
                    doc = procedure.doc_list
                    arguments = get_wrapper_arguments(procedure)
                    c_module += C_Wrapper(name, procedure.name, doc, arguments, module.name, getattr(procedure, "retvar", None))

            if len(c_module) > 0:
                c_modules.append(c_module)

    return c_modules


def get_wrapper_arguments(procedure: Procedure) -> C_Arguments:
    arguments = []
    for arg in procedure.args:
        dimension, extra_dim_args = get_dimension(arg)
        arguments.append(C_Argument(
            name=arg.name,
            doc=arg.doc_list,
            type=get_c_type(arg),
            dimension=dimension,
            intent=Intent(arg.intent),
            default_value=get_default_value(arg)
        ))
        arguments.extend(extra_dim_args)

    if (retvar := getattr(procedure, "retvar", None)) is not None:
        dimension, extra_dim_args = get_dimension(retvar)
        result_argument = C_Argument(
            name=retvar.name,
            doc=retvar.doc_list,
            type=get_c_type(retvar),
            dimension=dimension,
            intent=Intent.OUT,
            only_c_arg=True
        )
        arguments.append(result_argument)
        arguments.extend(extra_dim_args)

    return tuple(arguments)


def get_default_value(arg: FortranVariable) -> int | str | bool | float | None:
    # TODO: parse default value from doc_list[-1]
    if arg.optional:
        ...


def get_dimension(arg: FortranVariable) -> Tuple[Dimension, Arguments]:
    dims = arg.dimension
    for attrib in arg.attribs:
        if (match := re.match(r".*dimension\s*(?P<dimension>\([^\)]+\)).*", attrib)) is not None:
            dims = match.group("dimension")

    dims = dims[1:-1]
    dimension = []
    extra_dim_args = []

    # Handle length argument for characters
    if arg.full_type.startswith("character"):
        len_match = re.match(r".*\blen\s*=\s*(?:(?P<assumed_len>\*)|(?P<explicit_len>[^\),])).*", arg.full_type, re.IGNORECASE)
        if len_match is None:
            raise SyntaxError("Dummy arguments of type 'character' need to be declared with a named 'len' type paramter, like 'character(len=42)'")

        # if assumed length: add new length argument and add it to first extent of dimension
        if (strlen := len_match.group("assumed_len")) is not None:
            arg = C_Argument(
                name=f"{arg.name}_strlen",
                doc=[f" String length of `{arg.name}`"],
                type=C_Type.C_INT,
                dimension=(),
                intent=Intent.IN,
                only_c_arg=True
            )
            dimension.insert(0, arg.name)
            extra_dim_args.append(arg)

        # if explicit length: only set first extent
        else:
            dimension.insert(0, len_match.group("explicit_len"))

    if dims:
        # convert dims to dimension tuple
        for dim_idx, dim in enumerate(dims.split(","), 1):
            dim = dim.strip()

            # Usually not the case, but exists and may be more common in future, assumed-shape
            # create new explicit size argument for that (needs update if we decide to use 'ISO_Fortran_binding.h')
            if dim == ":":
                arg = C_Argument(
                    name=f"n_{arg.name}_elements_dim_{dim_idx}",
                    doc=[f" Size of the {dim_idx}. dimension/extent of `{arg.name}`"],
                    type=C_Type.C_INT,
                    dimension=(),
                    intent=Intent.IN,
                    only_c_arg=True
                )
                dimension.append(arg.name)
                extra_dim_args.append(arg)

            # if not assumed-shape, it is either assumed-size ('*') or explicit size, either constant or a variable
            # In both cases it is interoperable with our wrappers
            else:
                dimension.append(dim)

    return tuple(dimension), tuple(extra_dim_args)


def get_c_type(arg: FortranVariable) -> str:
    if arg.full_type.startswith("complex"):
        type_name = f"complex({arg.kind})"
    elif arg.full_type.startswith("character"):
        type_name = "character"
    else:
        type_name = arg.kind
    return C_Type(type_name)


def get_wrapper_name(procedure: Procedure, name_suffix: str) -> str:
    """
        Creates the name for the C wrapper of the passed procedure, respecting the expert/non-expert convention
    """

    match type(procedure).__name__:
        case "FortranSubroutine":
            # alloc routines don't get the alloc suffix in the C wrapper name
            if procedure.name.endswith("_alloc"):
                name = procedure.name.removesuffix("_alloc")
            # non-alloc routines get the _expert suffix if there is an alloc variant
            elif procedure.parent.find_child(f"{procedure.name}_alloc") is not None:
                name = f"{procedure.name}_expert"
            # if no alloc variant, the base name is simply the subroutine name
            else:
                name = procedure.name
        case "FortranFunction":
            name = procedure.name
        case "FortranModuleProcedureImplementation":
            name = procedure.name
        case _:
            raise TypeError("Procedure must be FortranFunction or FortranSubroutine")

    return name + name_suffix
