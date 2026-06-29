import re
from .utils import CodeGenerator, regex_escaped_preprocessor, warn, error, iter_parents
from enum import Enum

NAME_PAT = "(?i:[a-z][a-z_0-9]*)"
TYPE_PATS = {
    "component": "(?:procedure|proc|subroutine|function|interface|absinterface|block|type|file|module|submodule|program|namelist)",
    "item": "(?:absinterface|bound|common|constructor|final|function|interface|modproc|subroutine|type|variable)"
}
def ford_link_part_pat(name): return rf"(?P<{name}>{NAME_PAT})\s*(?:\(\s*(?P<{name}_type>{TYPE_PATS[name]})\s*\)\s*)?"


class FordLink(CodeGenerator):
    """class for Ford Links like [[mod(module):sub(subroutine)]]"""
    RE = re.compile(rf"(?P<match>\[\[\s*{ford_link_part_pat("component")}(?::\s*{ford_link_part_pat("item")})?\]\])")

    def __init__(self, parent: DocLine, component: str, component_type: str = None, item: str = None, item_type: str = None):
        self.component = component
        self.component_type = component_type
        self.item = item
        self.item_type = item_type
        self.parent = parent

    @property
    def name(self):
        return str(self)

    def __str__(self):
        component = self.component
        if self.component_type is not None:
            component += f"({self.component_type})"
        item = ""
        if self.item is not None:
            item = f":{self.item}"
            if self.item_type is not None:
                item += f"({self.item_type})"

        return f"[[{component}{item}]]"

    def copy(self):
        return FordLink(self.parent, self.component, self.component_type, self.item, self.item_type)


class TextAlignment(CodeGenerator, Enum):
    LEFT = 0
    RIGHT = 1
    CENTER = 2

    @classmethod
    def _missing_(cls, value):
        if type(value) is str:
            match = re.match(r"(?P<left>:)?-+(?P<right>:)?", value)

            if match is not None:
                match match.groups():
                    case (None, ":"):
                        return cls.RIGHT
                    case (":", ":"):
                        return cls.CENTER
                    case _:
                        return cls.LEFT

            raise SyntaxError(f"'{value}' is wrong markdown table alignment")

        raise KeyError(f"'{value}' invalid markdown table alignment")


class FordTable(CodeGenerator):
    """class for markdown tables in Ford comments"""
    TABLE_RE = re.compile(r"\| .* \|")
    ALIGN_RE = re.compile(r"\|:?-+:?(?:\|:?-+:?)*\|")

    MODE_ALIASES = ("Mode", "Method")
    MODE_VAR_RE = {
        alias: re.compile(rf"^{alias.upper()}_(?P<mode_name>[A-Z_0-9]+)$")
        for alias in MODE_ALIASES
    }

    def __init__(self, parent: DocList, header: Tuple[DocLine, ...], align: Tuple[TextAlignment, ...], rows: List[Tuple[DocLine, ...]]):
        if len(rows) == 0:
            raise SyntaxError(f"No rows for table")

        n_cols = len(header)
        if len(align) != n_cols or any(len(row) != n_cols for row in rows):
            raise SyntaxError(f"Mismatching column counts for table, perhaps missing padding for pipes? (Wrong: | x| y|, Right: | x | y |)")

        self.parent = parent
        self.n_cols = n_cols
        self.header = header
        self.align = align
        self.rows = rows
        self._handle_mode_table()

    @property
    def name(self):
        return f"ford doc table (doc line {self.header[0]})"

    def copy(self):
        return FordTable(self.parent, tuple(map(DocLine.copy, self.header)), self.align, list(tuple(map(DocLine.copy, row)) for row in self.rows))

    def _handle_mode_table(self):
        self.mode_vars = ()
        self.mode_alias = None

        # Check for mode argument -> Has header: | Method | Value |
        if self.n_cols == 2 and self.header[1].parts == ("Value",) and len(self.header[0].parts) == 1 and self.header[0].parts[0] in self.MODE_ALIASES:

            alias = self.header[0].parts[0]
            mode_var_re = self.MODE_VAR_RE[alias]

            mode_vars = []
            for _, mode_var in self.rows:
                parts = mode_var.parts
                # Value column must be Ford Link, pointing to the respective variable
                # Method column is just for description, so doesn't need handling
                if len(parts) == 1 and type(parts[0]) is FordLink:
                    link = parts[0]
                    mode_name = link.item
                    if mode_var_re.match(mode_name) is None:
                        return

                    mode_vars.append(link)
                # if not Ford Link, the Value column is allowed to be empty
                elif parts != (""):
                    return

            self.mode_alias = alias
            self.mode_vars = tuple(mode_vars)
            for mode_var in self.mode_vars:
                mode_var.mode_name = mode_var_re.match(mode_var.item).group("mode_name")

    @classmethod
    def mode_var_table_example(cls, alias):
        assert alias in cls.MODE_ALIASES, f"Expected '{alias}' to be in {cls.MODE_ALIASES}"
        desc = f"bla {alias.lower()}"
        value = f"[[bla_module(module):{alias.upper()}_BLA(variable)]]"
        return f"""
Expected markdown table in Ford comment, like:
!! | {alias.center(len(desc), " ")} | {"Value".center(len(value), " ")} |
!! |-{"-" * len(desc)}-|-{"-" * len(value)}-|
!! | {desc} | {value} |
!! ...
"""

    @classmethod
    def from_str_list(cls, parent: DocList, str_list: List[str, ...], line_idx: int):
        table = super(FordTable, cls).__new__(cls)

        line = str_list[line_idx].strip()

        if cls.TABLE_RE.match(line) is None:
            return None, line_idx

        line = line[1:-1]

        header = tuple(DocLine.from_str(table, l, line_idx) for l in line.split(" | "))
        n_cols = len(header)

        line_idx += 1
        line = str_list[line_idx].strip()
        if cls.ALIGN_RE.match(line) is None:
            error(SyntaxError, f"Expected markdown table align, like: |-----|-----| in line {line_idx + 1}", parent)
        line = line[1:-1]
        align = tuple(map(TextAlignment, line.split("|")))

        rows = []
        for line_idx in range(line_idx + 1, len(str_list)):
            line = str_list[line_idx].strip()
            if cls.TABLE_RE.match(line) is None:
                break
            line = line[1:-1]

            rows.append(tuple(DocLine.from_str(table, l, line_idx) for l in line.split(" | ")))

        cls.__init__(table, parent, header, align, rows)
        return table, min(len(str_list), line_idx + 1)


class DocLine:
    """class for lines in a DocList"""
    def __init__(self, parent: DocList, line_idx: int, parts: Tuple[str | FordLink | FordTable, ...]):
        self.parts = parts
        self.parent = parent
        self.idx = line_idx

    @property
    def name(self):
        return f"ford doc line {self.idx + 1}"

    @classmethod
    def from_str_list(cls, parent: DocList, str_list: List[str, ...], line_idx: int):
        line = str_list[line_idx].strip()

        return cls.from_str(parent, line, line_idx), line_idx + 1

    @classmethod
    def from_str(cls, parent: DocList, line: str, line_idx: int):
        doc_line = super(DocLine, cls).__new__(cls)

        line = line.strip()

        # TODO parse Ford Syntax (@note)
        parts = []
        last = 0
        for match in FordLink.RE.finditer(line):
            start, end = match.span()

            # text before the match
            if start > last:
                parts.append(line[last:start])

            # the match object itself
            parts.append(FordLink(doc_line, match.group("component"), match.group("component_type"), match.group("item"), match.group("item_type")))

            last = end

        # trailing text
        if last < len(line):
            parts.append(line[last:])

        cls.__init__(doc_line, parent, line_idx, tuple(parts))
        return doc_line

    def copy(self):
        return DocLine(self.parent, self.idx, tuple(part.copy() if type(part) is not str else part for part in self.parts))

    def __format__(self, spec):
        return "".join(format(p, spec) if type(p) is not str else p for p in self.parts)


DocListItem = str | DocLine | FordTable


class DocList(CodeGenerator):
    """Class for managing parsed Ford documentation. Each line is one element"""
    DM_RE = {
        "required_if_mode": re.compile(regex_escaped_preprocessor.expand(r"DM_REQUIRED_IF_MODE((?P<mode_var_name>.*), (?P<mode_var_module_name>.*), (?P<mode_var>.*))")),
        "result_size_is": re.compile(regex_escaped_preprocessor.expand(r"DM_RESULT_SIZE_IS((?P<n_results>.*))")),
        "output_from_auto": re.compile(regex_escaped_preprocessor.expand(r"DM_OUTPUT_FROM((?P<output_name>.*), (?P<proc_name>.*), (?P<mod_name>.*), AUTO)")),
        "output_from_just_info": re.compile(regex_escaped_preprocessor.expand(r"DM_OUTPUT_FROM((?P<output_name>.*), (?P<proc_name>.*), (?P<mod_name>.*), JUST_INFO)")),
        "default": re.compile(regex_escaped_preprocessor.expand(r"DM_DEFAULT((?P<default_val>.*))")),
    }

    SUPPORTED_TYPES = ("module", "procedure", "argument", "variable")

    def __init__(self, parent, doc_list: List[DocListItem, ...], unit_type: str):
        if unit_type not in self.SUPPORTED_TYPES:
            raise ValueError(f"type must be one of: {", ".join(self.SUPPORTED_TYPES)}")

        while len(doc_list) > 0 and doc_list[-1] == "":
            doc_list.pop()

        self.parent = parent
        self.name = "ford doc"
        self.doc_list = []
        line_idx = 0
        while line_idx < len(doc_list):
            if type(doc_list[line_idx]) is str:
                table, line_idx = FordTable.from_str_list(self, doc_list, line_idx)
                if table is not None:
                    line = table
                else:
                    line, line_idx = DocLine.from_str_list(self, doc_list, line_idx)
            else:
                line = doc_list[line_idx].copy()
                line_idx += 1

            self.doc_list.append(line)
        self.type = unit_type

        self.meta = {key: None for key in self.DM_RE}
        self.meta["optional_output"] = None
        self.meta["mode_table"] = None

        for i, line in enumerate(self.doc_list):
            if type(line) is FordTable and len(line.mode_vars) > 0:
                self.meta["mode_table"] = line
            for key, regex in self.DM_RE.items():
                if (match := regex.match(str(line))) is not None:
                    self.meta[key] = match

        # special syntax

    def copy(self):
        return DocList(self.parent, self.doc_list, self.type)

    @classmethod
    def from_fortran(cls, parent, unit: FortranModule | FortranSubroutine | FortranFunction | FortranVariable | FortranModuleProcedureImplementation):
        match type(unit).__name__:
            case "FortranModule":
                ty = "module"
            case "FortranSubroutine" | "FortranFunction":
                ty = "procedure"
            case "FortranVariable":
                ty = "argument"
            case _:
                raise TypeError(f"DocList doesn't support '{type(unit).__name__}'")

        doc_list = unit.doc_list
        if len(doc_list) == 0 or (len(doc_list) == 1 and doc_list[0] == ""):
            warn("No Ford documentation", unit)

        doc_list = [line.strip() for line in doc_list]

        return cls(parent, doc_list, ty)

    def __getitem__(self, idx):
        return self.doc_list[idx]

    def __setitem__(self, idx, value):
        if type(value) not in DocListItem.__args__:
            raise TypeError(f"Can only set DocList item to a value of one of the following types: {",".join(ty.__name__ for ty in DocListItem.__args__)}, not {type(other).__name__}")
        self.doc_list[idx] = value

    def __add__(self, other: str):
        if type(other) is not str:
            raise TypeError(f"Can only add str to DocList, not {type(other).__name__}")
        return DocList(self.parent, [*self.doc_list, DocLine.from_str(self, other, -1)], self.type)

    def __radd__(self, other: str):
        if type(other) is not str:
            raise TypeError(f"Can only add DocList to str, not {type(other).__name__}")
        return DocList(self.parent, [DocLine.from_str(self, other, -1), *self.doc_list], self.type)

    def __len__(self):
        return len(self.doc_list)
