from codegen.api.c_wrapper import C_Wrapper_Modules
from codegen.api.fortran import Modules
from codegen.c_wrapper import C_Wrapper_Serializer
from codegen.python import Python_Serializer
from codegen.rcpp import Rcpp_Serializer
import warnings
warnings.simplefilter("ignore")

c_interface_dir = "src/c_interface"

c_mods = C_Wrapper_Modules(Modules(exclude_directories=[c_interface_dir]))

C_Wrapper_Serializer.dump(c_mods, out_dir=c_interface_dir)
Python_Serializer.dump(c_mods, out_dir="python/tensor_omics")
Rcpp_Serializer.dump(c_mods, out_dir="rcpp/tensor_omics")
