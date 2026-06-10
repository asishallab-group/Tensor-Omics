from .api.utils import CodeGenerator, Serializer, Indentable
from .api.c_wrapper import C_Wrapper_Modules
from .api.fortran import Intent


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
                        string = "np.real64"
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
            case "doc":
                ndim = len(self.dimension)
                if ndim > 0:
                    if ndim == 1 and self.name == "character":
                        string = "str"
                    else:
                        string = f"ndarray[{format(self, "dtype")}] of shape {format(self.dimension, "shape")} in column-major layout (order='F')"
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

    def DocList(self, spec):
        return Indentable("\n".join(self.doc_list))

    def Dimension(self, spec):
        match spec:
            case "shape":
                string = f"({", ".join(self)})"
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
            case "doc":
                if self.type.intent is not Intent.OUT and self.is_dim_arg_for is None:
                    optional = ", optional" if self.optional else ""
                    string = f"""{self.name} : {format(self.type, "doc")}{optional}
{format(self.doc) >> INDENT}"""
            case "arglist":
                if self.type.intent is not Intent.OUT and self.is_dim_arg_for is None and not self.is_shape_arg:
                    string = self.name
                    if self.default_value is not None:
                        string += "=" + str(self.default_value)
            case "ensure_numpy_array":
                if self.is_dim_arg_for is None and len(self.type.dimension) > 0 and not self.is_shape_arg:
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
            case "allocate":
                if self.type.intent is Intent.OUT:
                    string = f"{self.name} = {format(self.type, "init_val")}"
                elif self.type.name == "character" and len(self.type.dimension) > 0:
                    string = f'{self.name} = {self.name}.astype({format(self.type, "dtype")}, order="F")'
                elif len(self.type.dimension) == 0 and self.is_dim_arg_for is None:
                    string = f"{self.name} = {format(self.type, "ctypes")}({self.name})"
            case "create_dim_args":
                parent = self.is_dim_arg_for
                if parent is not None:
                    if parent.type.name == "character" and parent.type.dimension[0] == self.name:
                        string = f"{self.name} = {parent.name}.dtype.itemsize // {parent.name}.dtype.alignment"
                    else:
                        string = f"{self.name} = {parent.name}.shape[{parent.type.dimension.index(self.name) - (parent.type.name == "character")}]"
            case "argtypes":
                string = format(self.type, "argtypes")
            case "call_args":
                if len(self.type.dimension) == 0:
                    if self.type.intent is not Intent.OUT:
                        string = f"ctypes.byref({format(self.type, "ctypes")}({self.name}))"
                    else:
                        string = f"ctypes.byref({self.name})"
                else:
                    string = self.name
            case "readonly":
                if self.type.intent is not Intent.IN:
                    if len(self.type.dimension) > 1 or (self.name != "character" and len(self.type.dimension) > 0):
                        string = f"{self.name}.setflags(write=False)"
        return Indentable(string)

    def C_Wrapper_Arguments(self, spec):
        match spec:
            case "doc" | "allocate" | "ensure_numpy_array" | "create_dim_args":
                string = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "argtypes" | "call_args" | "arglist":
                string = ",\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "readonly":
                string = "\n".join(filter(bool, (format(arg, spec) for arg in self)))
            case "return":
                outputs = tuple(arg for arg in self if arg.type.intent is not Intent.IN and arg.name != "ierr")
                n_outputs = len(outputs)
                string = "return "
                if n_outputs == 0:
                    string += "None"
                elif n_outputs == 1:
                    string += outputs[0].name
                else:
                    dict_elements = ",\n".join(f'"{arg.name}": {arg.name}' for arg in outputs)
                    string += f"""{{
{Indentable(dict_elements) >> INDENT}
}}"""

        return Indentable(string)

    def C_Wrapper(self, spec):
        return f'''def {self.orig_procedure.name}(
{format(self.arguments, "arglist") >> INDENT * 2}
        ):
    """
{format(self.doc) >> INDENT}

    Parameters
    ----------
{format(self.arguments, "doc") >> INDENT}
    """

    # ensure all array inputs are numpy arrays
{format(self.arguments, "ensure_numpy_array") >> INDENT}

    # extract dimension arguments
{format(self.arguments, "create_dim_args") >> INDENT}

    # Create temporaries and/or outputs
{format(self.arguments, "allocate") >> INDENT}

    # define ctypes interface
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

    @classmethod
    def dump(cls, c_wrapper_modules: C_Wrapper_Modules, out_file: str = "python/tensoromics_functions.py", lib_path: str = "build/libtensor-omics.so"):
        from pathlib import Path

        c_wrapper_modules.use(cls)
        out_file = Path(out_file)
        if not out_file.parent.is_dir():
            raise ValueError(f"out_file's directory '{out_file.parent}' does not exist")

        with open(out_file, "w") as tox:
            tox.write(f"""from error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("{lib_path}")
tox = ctypes.CDLL(dll_path)


{format(c_wrapper_modules)}""")
