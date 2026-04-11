"""
This script creates snippets for:
 - calling Fortran module routines in Fortran (for files in "src/")
 - defining interfacing Python/R functions that call the corresponding C wrappers from Fortran (for files in "python/" and "rcpp/")
 - calling the defined interfacing functions

The snippets are being written to "snippets/*_snippets.json"
"""

import glob
import re
import json
from sys import stderr
from typing import Set, Tuple, Dict, List, TypedDict, Literal

Body = List[str]
Description = str


class Snippet(TypedDict):
    prefix: str
    body: Body
    description: Description


SnippetKind = Literal["f42", "tox"]
Snippets = Dict[str, Snippet | SnippetKind]
InterfacingLanguage = Literal["Python", "R"]
FortranModuleName = str
FortranModules = Set[FortranModuleName]
ModuleName = str | Literal["f42:helper", "tox:helper"]
WrapperName = str
Wrappers = Set[WrapperName]
DefinitionKind = Literal["subroutine", "function"]
ArgList = List[str]
FuncDef = Tuple[str, ArgList]


def main():
    fortran_wrappers, fortran_modules = generate_fortran_snippets("src/**/*.[fF]90", ["src/config.F90", "src/safeguard.F90"])
    generate_interfacing_snippets("python/*.py", "Python", fortran_modules, fortran_wrappers={*fortran_wrappers})
    generate_interfacing_snippets("rcpp/*.R", "R", fortran_modules, fortran_wrappers={*fortran_wrappers})


def get_initial_snippet_dicts() -> Tuple[Dict[Literal["_kind"], SnippetKind], Dict[Literal["_kind"], SnippetKind]]:
    """
        Returns a tuple of two dictionaries. They have the "_kind" key referring to their snipped kind ('tox' or 'f42')
    """
    return {"_kind": "tox"}, {"_kind": "f42"}


def generate_fortran_snippets(file_pattern: str, ignored_files: List[str] = []) -> Tuple[Wrappers, FortranModules]:
    """
        Parses a fortran file and generates call snippets for each subroutine inside a module.

        Returns: Tuple[
            Wrappers: set of c wrappers outside the modules
            FortranModules: set of module names
        ]
    """
    tox_snippets, f42_snippets = get_initial_snippet_dicts()
    fortran_wrappers: Wrappers = set()
    fortran_modules: FortranModules = set()

    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                multiline: bool | None = None            # will tell if a subroutine/function definition takes multiple lines (because of too many args)
                module: str | None = None                # will hold the module name
                module_snippets: Snippets | None = None  # will hold the snippet dictionary for current module's kind (f42/tox)
                func_name: str | None = None             # will hold the name of the current parsed function/subroutine definition
                args: str | None = None                  # will hold the parsed arguments of the current parsed function/subroutine definition
                def_kind: DefinitionKind | None = None   # will hold the kind of a definition, either subroutine or function
                snippet_count: int = 0                   # will count the number of generated snippets
                is_wrapper_section: bool = False         # tells if parser is below module definition

                for line in file:
                    # Fortran is not case sensitive, so we neither
                    line: str = line.lower()

                    if re.search(r"^end (sub)?module\b", line.lower()) is not None:
                        is_wrapper_section = True

                    # before doing anything, find the module definition
                    if module is None:
                        module_def = re.search(r"^(sub)?module\b\s*\(?(?P<module_name>[a-z_0-9]+)", line)
                        if module_def is not None:
                            module = module_def.group("module_name")
                            fortran_modules.add(module.lower())
                            module_snippets = get_module_snippets(module, f42_snippets, tox_snippets)
                    else:
                        # if there wasn't a multiline definition before, check current line for a new definition
                        if not multiline:
                            func_def = re.search(r"\b(?P<def_kind>function|subroutine)\s+(?P<func_name>[a-z_0-9]\w*)\s*\((?P<args>[a-z_0-9,\s]*)(?P<continuation_char>&|\))?", line)
                            if func_def is not None:
                                def_kind, func_name, args, continuation_char = func_def.groups()
                                multiline = continuation_char == "&"

                        # current line continues a started subroutine/function definition
                        else:
                            more_args, continuation_char = re.search(r"\s*(?P<more_args>[a-z_0-9, ]*)(?P<continuation_char>&|\))?", line).groups()
                            args += more_args
                            multiline = continuation_char == "&"

                        # once a definition is done
                        if not multiline and multiline is not None:
                            # add interfacing subroutines (C wrappers) to the wrapper set
                            if is_wrapper_section:
                                if func_name[-1].lower() == "c":
                                    fortran_wrappers.add(func_name.lower())

                            # build call snippet for current definition
                            else:
                                # subroutines are called with `call foo(${1:arg1})`
                                if def_kind == "subroutine":
                                    call_str = "call"
                                    arg_list_string = arg_list_to_string(map(str.strip, args.split(",")))

                                # functions are called with `${1:result} = foo(${2:arg1})`
                                else:
                                    call_str = "${{1:result}} ="
                                    arg_list_string = arg_list_to_string(map(str.strip, args.split(",")), 2)

                                # add snippet
                                module_snippets[f"Fortran Call of {def_kind} '{func_name}'"] = {
                                    "prefix": f"{module_snippets["_kind"]}:{func_name}",
                                    "body": f"{call_str} {func_name}({arg_list_string})",
                                    "description": f"Call of Fortran {def_kind} '{func_name}' from module '{module}'"
                                }

                                # reset variables
                                func_name = None
                                args = None
                                def_kind = None
                                multiline = None

                                # increase snippet count
                                snippet_count += 1

                if module is None:
                    raise RuntimeError(f"No module found in file: {file_name}")
                if snippet_count == 0:
                    raise RuntimeError(f"No snippets generated for file '{file_name}'. If this is unexpected, either the syntax is wrong or there is overusage of multiline characters (&). If it is expected, add the file to the ignored files in the main function.")

    write_snippets(f42_snippets, tox_snippets, "Fortran")
    return fortran_wrappers, fortran_modules


def get_module_snippets(module_name: str, f42_snippets: Snippets, tox_snippets: Snippets) -> Snippets:
    """
        Module names need to start with either "f42_" or "tox_". Depending on the prefix, the respective snippets object is returned.

        Note that for helpers the prefix is "f42:" or "tox:". This is allowed as well.
    """
    if re.search(r'^f42[_:]', module_name):
        return f42_snippets
    if re.search(r'^tox[_:]', module_name):
        return tox_snippets
    else:
        raise RuntimeError(f"module name '{module_name}' does not start with 'tox_' or 'f42_'")


def generate_interfacing_snippets(file_pattern: str, lang: InterfacingLanguage, fortran_modules: FortranModules, ignored_files: List[str] = [], fortran_wrappers: Wrappers = set()) -> None:
    """
        Parses a python and R files and generates snippets for those lines having a preceding header.
        For function definitions, a call snippet is generated as well.
        For more information about headers, see `get_interfacing_metadata`.

        The snippets will be written to "snippets/${lang}_(f42|tox)_snippets.json"
    """

    tox_snippets, f42_snippets = get_initial_snippet_dicts()
    for file_name in glob.glob(file_pattern, recursive=True):
        if all(map(lambda f: f != file_name, ignored_files)):
            with open(file_name) as file:
                body: Body = []                                  # will hold the lines for the body of the snippet
                description: Description | None = None           # will hold the description of the snippet
                wrapped_func_name: WrapperName | None = None     # will hold the name of the Python/R function (or the name for a helper)
                module: ModuleName | None = None                 # will hold the module name parsed from snippet header
                module_snippets: Snippets | None = None          # either tox or f42 snippet dictionary, depending on handled snippet
                defined_function: bool | None = None             # tells if a function definition appeared below the snippet header to avoid multiple definitions in the same snippet
                is_void: bool = True                             # tells if the defined function returns something (if it has a return statement in it)

                for line_number, line in enumerate(file, 1):
                    # snippet header found
                    if line.startswith("#>"):
                        if line.startswith("#>skip snippets"):
                            break
                        update_interfacing_snippets(module_snippets, body, module, wrapped_func_name, description, is_void, defined_function, lang, fortran_wrappers, file_name, line_number)
                        is_void = True  # assume a void function
                        body = []
                        module, wrapped_func_name, description = get_interfacing_metadata(line, file_name, line_number)
                        defined_function = False

                        # guarantee that the specified module and subroutine in the snippet header really exist in Fortran (TODO: check if the subroutine is really related to the module)
                        if not module.endswith(":helper"):
                            if module.lower() not in fortran_modules:
                                raise error_in_line(f"Unknown module '{module}'", file_name, line_number)
                            if wrapped_func_name.lower() not in fortran_wrappers:
                                raise error_in_line(f"Unknown interfacing subroutine '{wrapped_func_name}' for module '{module}'", file_name, line_number, "(Could also be duplicate)")

                        # for the `.*helper-<name>` case where no function is defined in the snippet, treat it as if there was one, so no new one is defined
                        elif wrapped_func_name is not None:
                            defined_function = True

                        module_snippets = get_module_snippets(module, f42_snippets, tox_snippets)

                    # module is only None if there hasn't been a snippet header yet,
                    # so if there was a header before, extend the body with the current line
                    elif module is not None:
                        # check for a function definition
                        if re.search(r'(<-\s*function|^def)', line) is not None:
                            if defined_function:
                                error_in_line("Declaring a new function without snippet header", file_name, line_number)
                            defined_function = True

                        # treat line as body (without line break)
                        body.append(line[:-1])

                        # if function is still assumed as void, check for return statement
                        if is_void:
                            # check if return statement occurs
                            if (match := re.search(r"\breturn\b\s*(?P<return_value>\S+)?", line)) is not None:
                                return_value = match.group("return_value")

                                # if there follows text after the return statement, check if it is still some void stuff or not
                                if return_value is not None:
                                    is_comment = return_value.startswith("#")
                                    is_R_void = lang == "R" and re.search(r'^\(\s*invisible\s\(\sNULL\s\)\)', return_value) is not None
                                    if not (is_comment or is_R_void):
                                        is_void = False

                update_interfacing_snippets(module_snippets, body, module, wrapped_func_name, description, is_void, defined_function, lang, fortran_wrappers, file_name, line_number)

    # guarantee that each Fortran C wrapper has a related python snippet/function
    # TODO: once Rcpp is complete, remove the Python condition everywhere in this script
    if lang == "Python" and len(fortran_wrappers) > 0:
        raise RuntimeError(f"Missing Python wrappers for: {", ".join(fortran_wrappers)}")

    write_snippets(f42_snippets, tox_snippets, lang)


def write_snippets(f42_snippets: Snippets, tox_snippets: Snippets, lang: InterfacingLanguage):
    """
        Serializes snippets to file "snippets/${lang}_(f42|tox)_snippets.json"
    """
    for snippets in ["tox_snippets", "f42_snippets"]:
        with open(f"snippets/{lang}_{snippets}.json", "w") as f:
            snippets_dict = locals()[snippets]
            kind = snippets_dict["_kind"]
            del snippets_dict["_kind"]  # don't serialize the snippet kind
            json.dump(snippets_dict, f, indent=2)
            snippets_dict["_kind"] = kind  # re-add the snippet kind


def error_in_line(msg: str, file_name: str, line_number: int, post_msg: str | None = None):
    """
        Helper for raising errors consistently that give info about file name and line number
        msg-pattern: "${msg} in file '${filename}', line ${line_number}[: ${post_msg}]"
    """
    err_msg = f"{msg} in file '{file_name}', line {line_number}"
    if post_msg is not None:
        err_msg = f"{err_msg}: {post_msg}"
    raise RuntimeError(err_msg)


def get_interfacing_metadata(initial_line: str, file_name: str, line_number: int) -> Tuple[ModuleName, WrapperName | None, Description]:
    """
        Extract metadata from snippet headers (stored somewhere in initial_line).
        Supported Snippet headers:
         - "#> <fortran_module_name>:<fortran_c_wrapper_name>:<description>"
         - "#> (f42|tox)_helper(-<helper_name>)?:<description>"

        Pattern for names: [a-zA-Z_0-9]+

        Returns: Tuple[
            ModuleName
            WrapperName | None : only None for a helper definition without helper_name
            Description
        ]
    """
    metadata = re.search(r"^#> *(?P<module>[a-z_0-9A-Z]+):(?P<wrapped_func_name>[a-z_0-9A-Z]+):\s*(?P<description>.*)\s*", initial_line)
    if metadata is not None:
        return metadata.group("module"), metadata.group("wrapped_func_name"), metadata.group("description")

    metadata = re.search(r"^#> (?P<kind>f42|tox)_helper(-(?P<name>[a-zA-Z0-9_]+))?:\s*(?P<description>.*)\s*", initial_line)
    if metadata is not None:
        return f"{metadata.group("kind")}:helper", metadata.group("name"), metadata.group("description")

    error_in_line("Invalid snippet header", file_name, line_number, "Should match pattern: '#> <fortran_module_name>:<fortran_subroutine_name>:<description>' OR '#> (tox|f42)_helper[-<name>]:<description>' ")


def update_interfacing_snippets(snippets: Snippets, body: Body, module: ModuleName | None, wrapped_func_name: WrapperName, description: Description, is_void: bool, defined_function: str, lang: InterfacingLanguage, fortran_wrappers: Wrappers, file_name: str, line_number: int) -> None:
    """
        Generate snippets from definition metadata and add them to snippets dictionary
    """

    # module is None if there hasn't been a snippet header yet
    if module is not None:
        # if there was no function definition -> error
        if not defined_function:
            error_in_line("Could not inherit a name for the snippet, try '#> (f42|tox)_helper-<name>:<description>' syntax", file_name, line_number - len(body) - 1)

        # for a named helper, generate its snippet
        if module.endswith(":helper") and wrapped_func_name is not None:
            remove_unnecessary_lines(body)
            snippets[f"{lang} Helper to {description}"] = {
                "prefix": f"{snippets["_kind"]}:{lang[:2].lower()}_helper_{wrapped_func_name}",
                "body": body,
                "description": description
            }

        # generate snippets for function definitions
        else:
            func_name, args = get_interfacing_function_definition(body, lang)

            # helpers don't have a fortran wrapper, so their function name will be taken for the snippet
            if module.endswith(":helper"):
                wrapped_func_name = func_name
            # remove fortran wrapper from set to avoid duplicates (TODO: once Rcpp migration is done, remove Python condition)
            elif lang == "Python":
                fortran_wrappers.discard(wrapped_func_name.lower())

            snippets.update(generate_definition_snippet(body, description, module, wrapped_func_name, lang, snippets["_kind"]))
            snippets.update(generate_interfacing_call_snippet(func_name, module, wrapped_func_name, args, is_void, lang, snippets["_kind"]))


def remove_unnecessary_lines(body: Body) -> None:
    """
        Helper that removes all bottom lines of the body that are either empty or comments
    """

    inside_multiline_comment = False
    # remove last unnecessary lines, either empty (or spaces) or comments
    while (match := re.search(r'^\s*(#|(?P<ml_comment_open>""".*(?P<ml_comment_close>""")?)|$)', (last_line := body.pop()))) is not None or inside_multiline_comment:
        if match is not None and match.group("ml_comment_open") is not None:
            if match.group("ml_comment_close") is None:
                inside_multiline_comment = not inside_multiline_comment

    body.append(last_line)
    body.append("")


def generate_definition_snippet(body: Body, description: Description, module: ModuleName, wrapped_func_name: WrapperName, lang: InterfacingLanguage, snippet_kind: SnippetKind) -> Snippets:
    """
        Generate snippet for creation of a specific Python/R definition
    """

    if module.endswith(":helper"):
        key = f"{lang} Helper function '{wrapped_func_name}'"
        prefix = f"{snippet_kind}:{lang[:2].lower()}_helper_{wrapped_func_name}"
    else:
        key = f"{lang} setup for '{wrapped_func_name}' from module '{module}'"
        prefix = f"{snippet_kind}:{lang[:2].lower()}_setup_{wrapped_func_name}"

    remove_unnecessary_lines(body)

    return {
        key: {
            "prefix": prefix,
            "body": body,
            "description": description
        }
    }


def arg_list_to_string(args: ArgList, tab_index=1) -> str:
    """
        Helper that joins the args with a comma after wrapping the argument in a snippet tab stop "${i:<arg>}"
    """
    return ", ".join(f"${{{i}:{arg}}}" for i, arg in enumerate(args, tab_index))


def generate_interfacing_call_snippet(func_name: str, module: ModuleName, wrapped_func_name: WrapperName, args: ArgList, is_void: bool, lang: InterfacingLanguage, snippet_kind: SnippetKind) -> Snippets:
    """
        Generate call snippet for a Python/R definition
    """

    match lang:
        case "Python":
            assignment_operator = "="
        case "R":
            assignment_operator = "<-"
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")

    if module.endswith(":helper"):
        description: Description = f"Call of {lang} helper '{wrapped_func_name}'"
    else:
        description: Description = f"Call of {lang} wrapper for '{wrapped_func_name}' from module '{module}'"

    if not is_void:
        body: Body = [f"${{1:result}} {assignment_operator} "]
        arg_list_string = arg_list_to_string(args, 2)
    else:
        body: Body = [""]
        arg_list_string = arg_list_to_string(args, 1)

    body[0] += f"{func_name}({arg_list_string})"

    return {
        f"{lang} call of '{func_name}'": {
            "prefix": f"{snippet_kind}:{lang[:2].lower()}_{func_name}",
            "body": body,
            "description": description
        }
    }


def get_interfacing_function_definition(body: Body, lang: InterfacingLanguage) -> FuncDef:
    """
    Extract (function_name, [args]) from a snippet body.
    `body` is a list of lines (strings).
    """

    match lang:
        case "Python":
            return _get_python_definition(body)
        case "R":
            return _get_r_definition(body)
        case _:
            raise ValueError(f"Unsupported interfacing language: {lang}")


def _get_python_definition(body: Body) -> FuncDef:
    """
        Parses the body for a function definition and returns the function name and its arguments
    """

    # Remove inline comments and join lines
    cleaned: Body = []
    for line in body:
        # strip inline comments
        line = re.sub(r"#.*", "", line)
        cleaned.append(line)

    text = "".join(cleaned)

    # def foo(a,b,c):
    m = re.search(r"def\s+(?P<func_name>[A-Za-z_]\w*)\s*\((?P<args>.*?)\)\s*(->\s*\w*\s*)?:", text, re.S)
    if not m:
        raise ValueError("No Python function definition found in:\n\n", text)

    func_name: str = m.group("func_name")
    arglist: str = m.group("args").strip()

    args: ArgList = []
    if arglist:
        for a in arglist.split(","):
            a = a.strip()
            a = a.split("=")[0].strip()  # remove defaults
            a = a.lstrip("*")            # remove *args / **kwargs
            if a:
                args.append(a)

    return func_name, args


def _get_r_definition(body: Body) -> FuncDef:
    """
        Parses the body for a function definition and returns the function name and its arguments
    """

    # Remove inline comments and join lines
    cleaned: Body = []
    for line in body:
        line = re.sub(r"#.*", "", line)
        cleaned.append(line)

    text = "".join(cleaned)

    # foo <- function(a,b,c)
    m = re.search(r"(?P<func_name>[A-Za-z_]\w*)\s*<-\s*function\s*\((?P<args>.*?)\)", text, re.S)
    if not m:
        raise ValueError("No R function definition found")

    func_name: str = m.group("func_name")
    arglist: str = m.group("args").strip()

    args: ArgList = []
    if arglist:
        for a in arglist.split(","):
            a = a.strip()
            a = a.split("=")[0].strip()
            if a:
                args.append(a)

    return func_name, args


if __name__ == '__main__':
    main()
