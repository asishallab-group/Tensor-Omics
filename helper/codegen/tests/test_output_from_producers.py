"""`DM_OUTPUT_FROM` producers, and what calling one takes.

A producer input that is a shape the consumer's bindings read off its array.
"""

from codegen.diagnostics import DiagnosticBag
from codegen.ir.directives import Directives, OutputFrom, OutputFromMode
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent

from builders import integer, module, procedure, project, real


def auto(output, producer, producer_module):
    return Directives(output_from=OutputFrom(output, producer, producer_module, OutputFromMode.AUTO))


class TestAShapeTheBindingsReadOffAnArray:
    """A producer input that is a shape the consumer's bindings read off its array.

    R spelled it by name, `arr_shape`, which is not a variable there -- the shim reads it off
    `arr` -- so the call ended in "object 'arr_shape' not found". Python has the variable.
    """

    def test_r_reads_it_off_the_array_as_the_shim_does(self):
        from codegen.abi.c_abi import build_project
        from codegen.emit.r_wrapper import RWrapperEmitter

        counter = procedure(
            "count_vectors",
            integer("n_dims", doc="rank"),
            integer("arr_shape", Intent.IN, "(n_dims)", doc="a shape"),
            integer("n_vectors", Intent.OUT, doc="how many vectors"),
            integer("ierr", Intent.OUT, doc="error"),
        )
        norm = procedure(
            "norm",
            integer("n_elements", doc="elements"),
            real("arr", Intent.IN, "(n_elements)", doc="the array"),
            integer("n_dims", doc="rank"),
            integer("arr_shape", Intent.IN, "(n_dims)", doc="its shape"),
            integer("n_vectors", doc="vectors",
                    directives=auto("n_vectors", "count_vectors", "tox_hand")),
            integer("status", Intent.OUT, "(n_vectors)", doc="per vector"),
            integer("ierr", Intent.OUT, doc="error"),
        )
        bag = DiagnosticBag()
        source = project(module("tox_hand", counter, norm))
        analyse_project(source, bag)
        binding = build_project(source, bag)
        assert bag.errors == (), bag.render()
        wrapper = next(w for m in binding for w in m if w.stripped_name == "norm")

        text = RWrapperEmitter().function(wrapper)

        assert "count_vectors(arr_shape = if (is.null(dim(arr))) length(arr) else dim(arr))" in text
