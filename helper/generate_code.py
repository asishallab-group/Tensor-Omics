from codegen.c_wrappers import C_Modules, Project

C_Modules(Project()).dump("src")
