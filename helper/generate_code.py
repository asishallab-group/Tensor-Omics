from codegen.api.c_wrapper import C_Modules, Project
from codegen.c_wrapper import C_Wrapper_Serializer

C_Modules.use(C_Wrapper_Serializer)
c_mods = C_Modules(Project())
c_mods.dump("src")
