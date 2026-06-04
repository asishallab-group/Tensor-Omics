from codegen.c_wrappers import generate_c_module_code, collect_c_modules
from codegen.utils import get_project

generate_c_module_code(collect_c_modules(get_project()), "src")
