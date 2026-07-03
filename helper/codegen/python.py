from .api.utils import CodeGenerator, Serializer, Indentable
from .api.c_wrapper import C_Wrapper_Modules
from .api.fortran import Intent, Dimension


INDENT = 4


class Python_Serializer(Serializer):
    def Fortran_Type(self, spec):
        match spec:
            case "ctypes":
                match self.name:
                    case "integer":
                        string = "ctypes.c_int"
                    case "logical":
                        string = "ctypes.c_int"
                    case "real":
                        string = "ctypes.c_double"
                    case "complex":
                        string = "ctypes.c_double_complex"
                    case "character":
                        string = "ctypes.c_char"
            case "dtype":
                match self.name:
                    case "integer":
                        string = "np.int32"
                    case "logical":
                        string = "np.int32"
                    case "real":
                        string = "np.float64"
                    case "complex":
                        string = "np.complex128"
                    case "character":
                        string = f'f"S{{{format(self.dimension, "char_len")}}}"'
            case "argtypes":
                ndim = len(self.dimension)
                if ndim == 0:
                    string = f"ctypes.POINTER({format(self, "ctypes")})"
                else:
                    dtype = format(self, "dtype")
                    dtype = "" if not dtype else ", dtype=" + dtype
                    flags = "F_CONTIGUOUS" if ndim > 1 else "C_CONTIGUOUS"
                    if self.name == "character":
                        ndim -= 1
                    string = f"np.ctypeslib.ndpointer(ndim={ndim}, flags='{flags}'{dtype})"
            case "doc_params":
                ndim = len(self.dimension)
                if ndim > 0:
                    if ndim == 1 and self.name == "character":
                        string = "str"
                    else:
                        dimension = self.dimension[1:] if self.name == "character" else self.dimension
                        string = f"np.ndarray[{format(self, "dtype")}] of shape {format(dimension, "shape")} in column-major layout (order='F')"
                else:
                    match self.name:
                        case "integer":
                            string = "int"
                        case "complex":
                            string = "complex"
                        case "character":
                            string = "str"
                        case "logical":
                            string = "bool"
                        case "real":
                            string = "float"
            case "init_val":
                ndim = len(self.dimension)
                if ndim > 0:
                    if self.name == "character":
                        string = f"np.zeros({format(self.dimension[1:], "shape")}, dtype={format(self, "dtype")}, order='F')"
                    else:
                        string = f"np.empty({format(self.dimension, "shape")}, dtype={format(self, "dtype")}, order='F')"
                else:
                    string = f"{format(self, "ctypes")}(0)"

        return Indentable(string)

    def DocList(self, sep="\n"):
        return Indentable(sep.join(format(line) for line in self.doc_list))

    def Dimension(self, spec):
        match spec:
            case "shape":
                string = f"({", ".join(self)}"
                if len(self) == 1:
                    string += ","
                string += ")"
            case "char_len":
                string = self[0]

        return Indentable(string)

    def Intent(self, spec):
        raise NotImplementedError()

    def Procedure(self, spec):
        raise NotImplementedError()

    def Procedure_Argument(self, spec):
        raise NotImplementedError()

    def Procedure_Arguments(self, spec):
        raise NotImplementedError()

    def C_Wrapper_Argument(self, spec):
        string = ""
        match spec:
            case "doc_params":
                if self.type.intent is not Intent.OUT and len(self.is_dim_arg_for) == 0:
                    optional = ", optional" if self.optional else ""
                    in_place = ", modified in-place" if self.type.intent is Intent.INOUT else ""
                    string = f"""{self.name} : {format(self.type, "doc_params")}{in_place}{optional}
{format(self.doc_list) >> INDENT}"""
            case "doc_returns":
                if self.type.intent is Intent.OUT:
                    orig_dim = self.type.dimension
                    if (match := self.doc_list.meta["result_size_is"]) is not None:
                        self.type.dimension = Dimension((*orig_dim[:-1], match.group("n_results")))
                    string = f"""{self.name} : {format(self.type, "doc_params")}
{format(self.doc_list) >> INDENT}"""
                    self.type.dimension = orig_dim
            case "arglist":
                if self.type.intent is not Intent.OUT and len(tuple(arg for arg in self.is_dim_arg_for if arg.type.intent is not Intent.OUT)) == 0 and not self.is_shape_arg and self.is_mask_count_arg_for is None:
                    string = self.name
                    if self.optional:
                        string += "=" + str(self.default_value)
            case "ensure_numpy_array":
                if len(self.type.dimension) > 0 and self.type.intent is not Intent.OUT:
                    if len(self.is_dim_arg_for) == 0 and not self.is_shape_arg:
                        if self.type.intent is Intent.INOUT:
                            string = f'assert type({self.name}) is np.ndarray and {self.name}.flags.f_contiguous and {self.name}.dtype == {format(self.type, "dtype")}, "\'{self.name}\' must be column-major numpy array (order=\'F\')"'
                        else:
                            if self.type.name == "character":
                                string = f"{self.name} = np.asarray({self.name})"
                            else:
                                dtype = format(self.type, "dtype")
                                if dtype:
                                    dtype = ", dtype=" + dtype
                                mem_layout = "fortran" if len(self.type.dimension) > 1 else "contiguous"
                                string = f"{self.name} = np.as{mem_layout}array({self.name}{dtype})"

                    shape_arg = self.shape_arg
                    if shape_arg is not None:
                        string += f"\n{shape_arg.name} = np.ascontiguousarray({self.name}.shape, dtype=np.int32)"

                    if self.optional and self.default_value is None:
                        string = f"""if {self.name} is not None:
{Indentable(string) >> INDENT}
"""

            case "allocate":
                if self.type.intent is Intent.OUT:
                    string = f"{self.name} = {format(self.type, "init_val")}"
                elif self.type.name == "character" and len(self.type.dimension) > 0:
                    string = f'{self.name} = {self.name}.astype({format(self.type, "dtype")}, order="F")'
                elif len(self.type.dimension) == 0 and len(self.is_dim_arg_for) == 0:
                    string = f"{self.name} = {format(self.type, "ctypes")}({self.name})"

            case "create_dim_args":
                parents = tuple(arg for arg in self.is_dim_arg_for if arg.type.intent is not Intent.OUT)
                if len(parents) > 0:
                    parent = parents[0]
                    if parent.type.name == "character" and parent.type.dimension[0] == self.name:
                        string = f"{self.name} = {parent.name}.dtype.itemsize // {parent.name}.dtype.alignment"
                    else:
                        string = f"{self.name} = {parent.name}.shape[{parent.type.dimension.index(self.name) - (parent.type.name == "character")}]"
            case "create_mask_count_args":
                parent = self.is_mask_count_arg_for
                if parent is not None:
                    string = f"{self.name} = {parent.name}.sum(axis=-1)"
            case "argtypes":
                string = format(self.type, "argtypes")
                if self.optional and self.default_value is None:
                    # allow null
                    string = f"nullable({string})"
            case "call_args":
                if len(self.type.dimension) == 0:
                    if self.type.intent is not Intent.OUT and len(self.is_dim_arg_for) > 0:
                        string = f"ctypes.byref({format(self.type, "ctypes")}({self.name}))"
                    else:
                        string = f"ctypes.byref({self.name})"
                else:
                    string = self.name
            case "readonly":
                if self.type.intent is not Intent.IN and not self.is_temporary:
                    if len(self.type.dimension) > 1 or (self.name != "character" and len(self.type.dimension) > 0):
                        string = f"{self.name}.setflags(write=False)"
            case "return":
                if len(self.type.dimension) == 0:
                    string = f"{self.name}.value"
                elif (match := self.doc_list.meta["result_size_is"]) is not None:
                    string = f"{self.name}[..., :{match.group("n_results")}.value].copy()"
                else:
                    string = f"{self.name}"

        return Indentable(string)

    def C_Wrapper_Arguments(self, spec):
        match spec:
            case "doc_params" | "allocate" | "ensure_numpy_array" | "create_dim_args" | "create_mask_count_args":
                string = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "argtypes" | "call_args" | "arglist":
                string = ",\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "readonly":
                string = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "return":
                outputs = tuple(arg for arg in self if arg.type.intent is Intent.OUT and arg.name != "ierr" and not arg.is_temporary)
                n_outputs = len(outputs)
                string = "return "
                if n_outputs == 0:
                    string += "None"
                elif n_outputs == 1:
                    string += format(outputs[0], "return")
                else:
                    dict_elements = ",\n".join(f'"{arg.name}": {format(arg, "return")}' for arg in outputs)
                    string += f"""{{
{Indentable(dict_elements) >> INDENT}
}}"""
            case "doc_returns":
                outputs = tuple(arg for arg in self if arg.type.intent is Intent.OUT and arg.name != "ierr")
                outputs_doc = tuple(format(arg, "doc_returns") for arg in outputs)
                n_outputs = len(outputs)
                if n_outputs == 0:
                    string = "None"
                elif n_outputs == 1:
                    string = outputs_doc[0]
                else:
                    dict_elements = ",\n".join(outputs_doc)
                    string = f"""results : dict
{Indentable(dict_elements) >> INDENT}
"""

        return Indentable(string)

    def C_Wrapper(self, spec):
        summary = self.orig_procedure.meta.summary
        if summary:
            summary = f"""\n    {summary}\n\n"""
        else:
            summary = ""
        if len(self.doc_list) > 0:
            notes = f"""
    Notes
    -----
{format(self.doc_list) >> INDENT}
"""
        else:
            notes = ""

        return f'''def {self.stripped_name}(
{format(self.arguments, "arglist") >> INDENT * 2}
        ):
    """{summary}
    Parameters
    ----------
{format(self.arguments, "doc_params") >> INDENT}

    Returns
    -------
{format(self.arguments, "doc_returns") >> INDENT}
{notes}    """

    # ensure all array inputs are numpy arrays
{format(self.arguments, "ensure_numpy_array") >> INDENT}

    # extract dimension arguments
{format(self.arguments, "create_dim_args") >> INDENT}
{format(self.arguments, "create_mask_count_args") >> INDENT}

    # Create temporaries and/or outputs
{format(self.arguments, "allocate") >> INDENT}

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {{'from_param': from_param}})

    tox.{self.name}.argtypes = (
{format(self.arguments, "argtypes") >> 2 * INDENT}
    )
    tox.{self.name}.restype = None

    tox.{self.name}(
{format(self.arguments, "call_args") >> 2 * INDENT}
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
{format(self.arguments, "readonly") >> INDENT}

{format(self.arguments, "return") >> INDENT}
'''

    def C_Wrapper_Module(self, spec):
        return "\n\n".join(format(c_wrapper) for c_wrapper in self.c_wrappers)

    def C_Wrapper_Modules(self, spec):
        return "\n\n".join(format(mod) for mod in self)

    def Error_Handling(self, spec):
        return f"""ERROR_CODES = {{
{Indentable(",\n".join(f'{code}: "{format(doc_list, " ").strip()}"' for code, doc_list in self.error_codes.values() if code != 0)) >> INDENT}
}}


def check_err_code(ierr: int):
    if ierr == 0:
        return

    msg = ERROR_CODES.get(ierr, f"Unmapped error code: {{ierr}}")
    raise RuntimeError(msg)
"""

    @classmethod
    def dump(cls, c_wrapper_modules: C_Wrapper_Modules, out_dir: str = "python/tensor_omics", lib_path: str = "build/libtensor-omics.so"):
        from pathlib import Path
        from os import mkdir

        c_wrapper_modules.use(cls)
        out_dir = Path(out_dir)
        if not out_dir.parent.is_dir():
            raise ValueError(f"out_dir's directory '{out_dir.parent}' does not exist")

        if out_dir.is_dir():
            from shutil import rmtree
            rmtree(out_dir)

        mkdir(out_dir)

        with open(out_dir.joinpath("error_handling.py"), "w") as error_handling:
            error_handling.write(format(c_wrapper_modules.error_handling))

        for module in c_wrapper_modules:
            with open(out_dir.joinpath(module.stripped_name + ".py"), "w") as module_file:
                module_file.write(f"""from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("{lib_path}")
tox = ctypes.CDLL(dll_path)


{format(module)}""")
