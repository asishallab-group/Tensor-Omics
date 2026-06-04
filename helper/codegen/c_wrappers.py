from .ford_api import Project, Procedure, ProcedureArgument, DocList, Dimension, Intent, Fortran_Type
from .utils import Indentable
from typing import List, Tuple
from enum import Enum
import re
from pathlib import Path
import os


NAME_SUFFIX = "_c"
INDENT = 4


class C_Argument:
    """Includes all relevant information about an argument, meaning its type, dimension, name, docstring and default value if optional"""
    def __init__(self, name: str, doc: DocList, type: Fortran_Type, dimension: Dimension, intent: Intent, default_value=None, only_c_arg=False):
        self.name = name
        self.doc = doc
        self.type = type
        self.dimension = dimension
        self.intent = intent
        self.default_value = default_value
        self.only_c_arg = only_c_arg

    @classmethod
    def create_dim_arg(cls, name: str, orig_name: str, doc: str, only_c_arg=True):
        return cls(
            name=name,
            doc=[f" {doc}"],
            type=Fortran_Type("integer", kind="int32"),
            dimension=Dimension(),
            intent=Intent.IN,
            only_c_arg=only_c_arg
        )

    @classmethod
    def from_proc_arg(cls, arg: ProcedureArgument) -> Tuple[Self, Tuple[C_Argument]]:
        dimension = list(arg.dimension)
        extra_dim_args = []

        for dim_idx, dim in enumerate(arg.dimension):
            dim_name_suffix = f"_dim_{dim_idx + 1}" if len(arg.dimension) > 1 else ""
            if dim == ":":
                extra_arg = cls.create_dim_arg(
                    f"n_{arg.name}_elements" + dim_name_suffix,
                    arg.name,
                    f" Size of the {dim_idx + 1}. dimension/extent of `{arg.name}`"
                )
                dimension[dim_idx] = extra_arg.name
                extra_dim_args.append(extra_arg)

        if arg.type.name == "character":
            extra_arg = cls.create_dim_arg(
                f"{arg.name}_strlen",
                arg.name,
                f"String length of '{arg.name}'"
            )
            extra_dim_args.append(extra_arg)
            dimension.insert(0, extra_arg.name)

        argument = cls(
            name=arg.name,
            doc=arg.doc_list,
            type=arg.type,
            dimension=Dimension(dimension),
            intent=arg.intent,
            default_value=arg.default_value
        )

        return argument, tuple(extra_dim_args)

    def __format__(self, spec):
        match spec:
            case "dummy":
                arg_str = f"{format(self.type, "C")}, {self.dimension}{self.intent}, target :: {self.name}"
                for line in self.doc:
                    arg_str += f"\n    !!{line}"
            case "name":
                arg_str = self.name
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")
        return Indentable(arg_str)


class C_Arguments(tuple):
    """Collection of C_Argument objects"""
    def __new__(cls, procedure: Procedure):
        arguments = []
        for proc_arg in procedure.args:
            arg, extra_dim_args = C_Argument.from_proc_arg(proc_arg)
            arguments = [*extra_dim_args, *arguments, arg]

        if procedure.retvar is not None:
            retvar = procedure.retvar
            result_argument = C_Argument(
                name=retvar.name,
                doc=DocList.from_fortran(retvar),
                type=Fortran_Type.from_fortran_variable(retvar),
                dimension=Dimension.from_fortran_variable(procedure.retvar),
                intent=Intent.OUT,
                only_c_arg=True
            )
            arguments = [*extra_dim_args, *arguments, result_argument]

        return super(cls, cls).__new__(cls, arguments)

    def __format__(self, spec):
        match spec:
            case "dummy":
                formatted = "\n".join(format(arg, "dummy") for arg in self)
            case "arglist":
                formatted = ", ".join(arg.name for arg in self)
            case "null_validation":
                formatted = "M_CHECK_IERR_NON_NULL"
                for arg in self:
                    if arg.name != "ierr":
                        formatted += f"\nM_CHECK_NON_NULL({arg.name})"
            case _:
                raise ValueError(f"Unsupported format spec '{spec}'")

        return Indentable(formatted)


class C_Wrapper:
    """Includes all relevant information about a C wrapper, meaning its C name, arguments and argument docstrings"""
    def __init__(self, procedure: Procedure, module_name: str):
        self.doc = procedure.doc_list
        self.module_name = module_name
        self.orig_proc_name = procedure.name
        self.orig_proc_type = procedure.type
        self.orig_proc_ford_link = f"[[{self.module_name}(module):{self.orig_proc_name}({self.orig_proc_type})]]"
        self.orig_proc_call = format(procedure, "call")

        match self.orig_proc_type:
            case "subroutine":
                # alloc routines don't get the alloc suffix in the C wrapper name
                if procedure.name.endswith("_alloc"):
                    name = procedure.name.removesuffix("_alloc")
                # non-alloc routines get the _expert suffix if there is an alloc variant
                elif procedure.parent.find_child(f"{procedure.name}_alloc") is not None:
                    name = f"{procedure.name}_expert"
                # if no alloc variant, the base name is simply the subroutine name
                else:
                    name = procedure.name
            case "function":
                name = procedure.name

        self.name = name + NAME_SUFFIX
        self.arguments = C_Arguments(procedure)

    def __str__(self):
        wrapper = f"""!> C-wrapper for {self.orig_proc_ford_link}
!|{"\n!|".join(self.doc)}
subroutine {self.name}({format(self.arguments, "arglist")}) bind(C, name="{self.name}")
    use {self.module_name}, only: {self.orig_proc_name}
{format(self.arguments, "dummy") >> INDENT}

{format(self.arguments, "null_validation") >> INDENT}

{self.orig_proc_call >> INDENT}
end subroutine {self.name}"""

        return Indentable(wrapper)


class C_Module:
    """Includes a bunch of wrappers that should appear in one c interfacing module"""
    def __init__(self, name: str, orig_mod_name: str, doc: DocList):
        self.name = name
        self.orig_mod_name = orig_mod_name
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

!> Module for C-wrappers for [[{self.orig_mod_name}(module)]]
!|{"\n!|".join(self.doc)}
module {self.name}
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

{"\n\n    ".join(str(c_wrapper) >> INDENT for c_wrapper in self)}

end module {self.name}"""


class C_Modules(tuple):
    """Collection of C_Module objects"""
    def __new__(cls, project: Project):
        """
            Creates and Collects all C wrappers for the subroutines/functions in the Project that have the 'category: C-interface' meta-tag and returns them.
            The output will only have non-empty modules.
        """
        c_modules = []

        for module in project.modules:
            if not module.name.endswith(NAME_SUFFIX):
                c_module = C_Module(module.name + NAME_SUFFIX, module.name, DocList.from_fortran(module))
                for procedure in map(Procedure, module.routines):
                    if procedure.meta.category == "C-interface":
                        c_module += C_Wrapper(procedure, module.name)

                if len(c_module) > 0:
                    c_modules.append(c_module)

        return super(cls, cls).__new__(cls, c_modules)

    def dump(self, out_dir: str):
        out_dir = Path(out_dir)
        if not out_dir.is_dir():
            raise ValueError("out_dir muste be a valid directory path")

        c_wrapper_dir = out_dir.joinpath("c_interface")
        if c_wrapper_dir.is_dir():
            from shutil import rmtree
            rmtree(c_wrapper_dir)

        os.mkdir(c_wrapper_dir)

        for module in self:
            module_file_path = c_wrapper_dir.joinpath(module.name + ".F90")

            with open(module_file_path, "w") as module_file:
                module_file.write(str(module))
