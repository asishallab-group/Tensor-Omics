from .fortran import Module, Procedure, Procedure_Argument, DocList, Dimension, Intent, Fortran_Type
from .utils import CodeGenerator
from typing import List, Tuple
from enum import Enum
import re
from pathlib import Path
import os


NAME_SUFFIX = "_c"


class C_Wrapper_Argument(CodeGenerator):
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
    def from_proc_arg(cls, arg: Procedure_Argument) -> Tuple[Self, Tuple[C_Wrapper_Argument]]:
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


class C_Wrapper_Arguments(CodeGenerator, tuple):
    """Collection of C_Wrapper_Argument objects"""
    def __new__(cls, procedure: Procedure):
        arguments = []
        for proc_arg in procedure.args:
            arg, extra_dim_args = C_Wrapper_Argument.from_proc_arg(proc_arg)
            arguments = [*extra_dim_args, *arguments, arg]

        if procedure.retvar is not None:
            retvar = procedure.retvar
            result_argument = C_Wrapper_Argument(
                name=retvar.name,
                doc=DocList.from_fortran(retvar),
                type=Fortran_Type.from_fortran_variable(retvar),
                dimension=Dimension.from_fortran_variable(procedure.retvar),
                intent=Intent.OUT,
                only_c_arg=True
            )
            arguments = [*extra_dim_args, *arguments, result_argument]

        return super(cls, cls).__new__(cls, arguments)


class C_Wrapper_Wrapper(CodeGenerator):
    """Includes all relevant information about a C wrapper, meaning its C name, arguments and argument docstrings"""
    def __init__(self, procedure: Procedure):
        self.doc = procedure.doc_list
        self.orig_procedure = procedure

        match procedure.type:
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
        self.arguments = C_Wrapper_Arguments(procedure)


class C_Wrapper_Module(CodeGenerator):
    """Includes a bunch of wrappers that should appear in one c interfacing module"""
    def __init__(self, module: Module):
        self.name = module.name + NAME_SUFFIX
        self.orig_module = module
        self.doc = module.doc_list
        self.c_wrappers = []
        for procedure in module.procedures:
            if procedure.meta.category == "C-interface":
                self.c_wrappers.append(C_Wrapper_Wrapper(procedure))

    def __iter__(self):
        return iter(self.c_wrappers)

    def __len__(self):
        return len(self.c_wrappers)


class C_Wrapper_Modules(CodeGenerator, tuple):
    """Collection of C_Wrapper_Module objects"""
    def __new__(cls, modules: Modules):
        """
            Creates and Collects all C wrappers for the subroutines/functions in the Project that have the 'category: C-interface' meta-tag and returns them.
            The output will only have non-empty modules.
        """
        c_modules = []

        for module in modules:
            if not module.name.endswith(NAME_SUFFIX):
                c_module = C_Wrapper_Module(module)
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
                module_file.write(format(module))
