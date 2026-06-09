from codegen.api.c_wrapper import C_Wrapper_Modules
from codegen.api.fortran import Modules
from codegen.c_wrapper import C_Wrapper_Serializer
from codegen.python import Python_Serializer

c_mods = C_Wrapper_Modules(Modules())
# C_Wrapper_Modules.use(C_Wrapper_Serializer)
# c_mods.dump("src")
c_mods.use(Python_Serializer)
Python_Serializer.dump(c_mods, out_file="python/tox.py")
