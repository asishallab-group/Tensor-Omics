from codegen.api.c_wrapper import C_Wrapper_Modules
from codegen.api.fortran import Modules
from codegen.c_wrapper import C_Wrapper_Serializer
from codegen.python import Python_Serializer

c_mods = C_Wrapper_Modules(Modules())

C_Wrapper_Serializer.dump(c_mods, out_dir="src")
Python_Serializer.dump(c_mods, out_dir="python/tensor_omics")
