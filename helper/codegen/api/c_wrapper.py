from .fortran import Module, Procedure, Procedure_Argument, DocList, Dimension, Intent, Fortran_Type, Procedure
from .utils import CodeGenerator, regex_escaped_preprocessor
from typing import List, Tuple
from enum import Enum
import re

NAME_SUFFIX = "_c"

MODE_TABLE_HEAD_RE = re.compile(r"\s*\|\s*(?:Mode|Method)\s*\|\s*Value\s*\|\s*")
MODE_TABLE_ALIGN_RE = re.compile(r"\s*\|:?-+:?\|:?-+:?\|\s*")
MODE_FORD_LINK_RE = re.compile(r"\[\[(?P<module_name>[a-z_]+)\(module\):(?P<mode_name>[A-Z_]+)\(variable\)\]\]")
MODE_TABLE_ROW_RE = re.compile(r"\s*\|[^\|]+\|\s*" + MODE_FORD_LINK_RE.pattern + r"\s*\|\s*")
MODE_VAR_PREFIX_RE = re.compile(r"^(?:MODE|METHOD)_")

Variable_Name = str
Mode_String = str
Module_Name = str


class C_Wrapper_Argument(CodeGenerator):
    """Includes all relevant information about an argument, meaning its type, dimension, name, docstring and default value if optional"""
    def __init__(self, name: str, doc: DocList, type: Fortran_Type, is_temporary: bool, c_wrapper: C_Wrapper, mode_vars: Tuple[Tuple[Variable_Name, Module_Name, Mode_String], ...] = None, default_value=None, optional=False, only_c_arg=False, is_dim_arg_for: Tuple[C_Wrapper_Argument | Procedure_Argument] = (), shape_arg: C_Wrapper_Argument | Procedure_Argument = None, is_shape_arg=False, is_mask_count_arg_for: Procedure_Argument | None = None):
        self.name = name
        self.doc_list = doc
        self.type = type
        self.default_value = default_value
        self.optional = optional
        self.only_c_arg = only_c_arg
        self.is_temporary = is_temporary
        self._is_dim_arg_for = is_dim_arg_for
        self._shape_arg = shape_arg
        self._is_mask_count_arg_for = is_mask_count_arg_for
        self.parent = c_wrapper
        self.is_shape_arg = is_shape_arg
        self.mode_vars = mode_vars

    @classmethod
    def create_dim_arg(cls, name: str, doc: str, c_wrapper: C_Wrapper, only_c_arg=True, is_dim_arg_for=None):
        return cls(
            name=name,
            doc=DocList([f" {doc}"], type="argument"),
            type=Fortran_Type("integer", intent=Intent.IN, kind="int32"),
            only_c_arg=only_c_arg,
            is_temporary=False,
            is_dim_arg_for=is_dim_arg_for,
            c_wrapper=c_wrapper
        )

    @classmethod
    def from_proc_arg(cls, arg: Procedure_Argument, c_wrapper: C_Wrapper) -> Tuple[Self, Tuple[C_Wrapper_Argument]]:
        dimension = list(arg.type.dimension)
        extra_dim_args = []

        shape_arg = arg.shape_arg

        if shape_arg is None:
            for dim_idx, dim in enumerate(arg.type.dimension):
                dim_name_suffix = f"_dim_{dim_idx + 1}" if len(arg.type.dimension) > 1 else ""
                if dim == ":":
                    extra_arg = cls.create_dim_arg(
                        f"n_{arg.name}_elements" + dim_name_suffix,
                        f"Size of the {dim_idx + 1}. dimension/extent of `{arg.name}`",
                        c_wrapper=c_wrapper
                    )
                    dimension[dim_idx] = extra_arg.name
                    extra_dim_args.append(extra_arg)
        else:
            # use assumed size for arguments whose shape is specified by a different argument
            dimension = ["*"]

        if arg.type.name == "character":
            if arg.type.len == "*":
                extra_arg = cls.create_dim_arg(
                    f"{arg.name}_strlen",
                    f"String length of '{arg.name}'",
                    c_wrapper=c_wrapper
                )
                extra_dim_args.insert(0, extra_arg)
                dimension.insert(0, extra_arg.name)
            elif arg.type.len == ":":
                raise AssertionError(f"Deferred length not allowed for '{arg.name}' in '{arg.parent.name}' of '{arg.parent.parent.name}'")
            else:
                dimension.insert(0, arg.type.len)

        doc_list = arg.doc_list

        if arg.is_mode_arg:
            state = "doc"
            converted_doc_list = []
            mode_vars = []
            for doc in doc_list:
                match state:
                    case "doc":
                        if MODE_TABLE_HEAD_RE.match(doc):
                            state = "mode_spec_align"
                    case "mode_spec_align":
                        if MODE_TABLE_ALIGN_RE.match(doc):
                            state = "mode_spec_rows"
                        else:
                            raise AssertionError(f"Expected some pattern like '|------|------|' for '{arg.name}' in '{arg.parent.name}' in '{arg.parent.parent.name}'")
                    case "mode_spec_rows":
                        if (match := MODE_TABLE_ROW_RE.match(doc)) is not None:
                            mode_name = match.group("mode_name")
                            module_name = match.group("module_name")
                            mode_string = MODE_VAR_PREFIX_RE.sub("", mode_name).lower()
                            mode_vars.append((mode_name, module_name, mode_string))
                            doc = MODE_FORD_LINK_RE.sub(f' "{mode_string}"  ', doc)
                        else:
                            state = "done"

                converted_doc_list.append(doc)

            assert len(mode_vars) > 0, f"""Found mode argument '{arg.name}' in '{arg.parent.name}' in '{arg.parent.parent.name}'. In Ford doc comment, expected markdown table like:
!! |    Mode   |   Value  |
!! |-----------|----------|
!! | bla mode  | MODE_BLA |
!! ...

Possible variants:
 - Mode->Method, MODE_BLA->METHOD_BLA
"""
            mode_vars = tuple(mode_vars)
            doc_list = DocList(converted_doc_list, "argument")
            assert arg.type.name == "integer", f"Found mode argument '{arg.name}' in '{arg.parent.name}' in '{arg.parent.parent.name}'. Must be integer, got: {arg.type.name}"
            max_mode_len = str(max(len(mode_str) for _, _, mode_str in mode_vars))
            arg_type = Fortran_Type("character", intent=Intent.IN, dimension=Dimension([max_mode_len, *dimension]), length=max_mode_len)
        else:
            arg_type = Fortran_Type(arg.type.name, intent=arg.type.intent, kind=arg.type.kind, dimension=Dimension(dimension), length=arg.type.len)
            mode_vars = None

        argument = cls(
            name=arg.name,
            doc=doc_list,
            type=arg_type,
            is_temporary=arg.is_temporary,
            is_dim_arg_for=arg.is_dim_arg_for,
            is_shape_arg=arg.is_shape_arg,
            shape_arg=shape_arg,
            default_value=arg.default_value,
            optional=arg.optional,
            c_wrapper=c_wrapper,
            mode_vars=mode_vars,
            is_mask_count_arg_for=arg.is_mask_count_arg_for
        )

        is_dim_arg_for = (argument,)
        for dim_arg in extra_dim_args:
            dim_arg._is_dim_arg_for = is_dim_arg_for

        return argument, tuple(extra_dim_args)

    @property
    def is_dim_arg_for(self):
        dim_args = list(self._is_dim_arg_for)
        for i_arg, arg in enumerate(dim_args):
            if type(arg) is Procedure_Argument:
                for c_arg in self.parent.arguments:
                    if c_arg.name == arg.name:
                        dim_args[i_arg] = arg

        return tuple(dim_args)

    @property
    def is_mask_count_arg_for(self):
        arg = self._is_mask_count_arg_for
        if type(arg) is Procedure_Argument:
            for c_arg in self.parent.arguments:
                if c_arg.name == arg.name:
                    return c_arg

    @property
    def shape_arg(self):
        arg = self._shape_arg
        if type(arg) is Procedure_Argument:
            for c_arg in self.parent.arguments:
                if c_arg.name == arg.name:
                    return c_arg


class C_Wrapper_Arguments(CodeGenerator, tuple):
    """Collection of C_Wrapper_Argument objects"""
    def __new__(cls, procedure: Procedure, c_wrapper: C_Wrapper):
        arguments = []
        ierr_seen = False
        for proc_arg in procedure.args:
            arg, extra_dim_args = C_Wrapper_Argument.from_proc_arg(proc_arg, c_wrapper)
            arguments.append(arg)
            arguments.extend(extra_dim_args)
            ierr_seen = ierr_seen or arg.name == "ierr"

        if procedure.retvar is not None:
            retvar = procedure.retvar
            ret_type = Fortran_Type.from_fortran_variable(retvar)
            ret_type.intent = Intent.OUT
            result_argument = C_Wrapper_Argument(
                name=retvar.name,
                doc=DocList.from_fortran(retvar),
                type=ret_type,
                only_c_arg=True,
                is_temporary=False,
                c_wrapper=c_wrapper
            )
            arguments.append(result_argument)

        if not ierr_seen:
            arguments.append(C_Wrapper_Argument(
                name="ierr",
                doc=DocList([" Error code"], "argument"),
                type=Fortran_Type("integer", Intent.OUT, kind="int32"),
                only_c_arg=True,
                is_temporary=False,
                c_wrapper=c_wrapper
            ))

        return super().__new__(cls, arguments)


class C_Wrapper(CodeGenerator):
    """Includes all relevant information about a C wrapper, meaning its C name, arguments and argument docstrings"""
    def __init__(self, procedure: Procedure):
        self.doc_list = procedure.doc_list
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
        self.arguments = C_Wrapper_Arguments(procedure, self)


class C_Wrapper_Module(CodeGenerator):
    """Includes a bunch of wrappers that should appear in one c interfacing module"""
    def __init__(self, module: Module):
        self.name = module.name + NAME_SUFFIX
        self.orig_module = module
        self.doc_list = module.doc_list
        self.c_wrappers = []
        for procedure in module.procedures:
            if procedure.meta.category == "C-interface":
                self.c_wrappers.append(C_Wrapper(procedure))

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

        mods = super().__new__(cls, c_modules)
        mods.project = modules.project
        return mods
