This directory includes frequently used or testable units of logic reused across development stages. 
- Snippets should be easy to create and use. The goal is to give the user access to the subroutine names along with their respective arguments, and nothing more. Example:
```
{
    "Call to subroutine_name": {
        "prefix": "tox|f42:subroutine_name",
        "body": [
        "! Brief explanation on what the subroutine does",
        "call subroutine_name(arg1, arg2, arg3, arg4)"
        ],
        "description": "Insert a call to subroutine_name with brief explanation"
    }
}
```
---

### How to Add Snippets in VSCode

To install a snippet in Visual Studio Code:

1. Open the command palette (`Ctrl+Shift+P` or `Cmd+Shift+P`).
2. Search for:  
   `Preferences: Configure Snippets`
3. Choose a global file (e.g. `global.code-snippets`) or the language-specific file for Fortran (e.g. `fortran.json`).
4. Paste the snippet block inside the JSON file.

---

### Snippets explanation

#### `f42:do-parallel`

Start typing `f42:do-parallel` and let the snippet autocomplete.

It includes tab stops for easy customization:
- `${1:i}` → loop index variable
- `${2:N}` → name of the loop limit variable
- `${3:10}` → initial value assigned to `N`
- Loop body → the cursor jumps there after filling in the loop setup

The snippet expands into a unified parallel loop structure that automatically selects the appropriate backend based on compiler flags:
- Coarray
- OpenMP
- Serial fallback

This allows backend-independent parallel execution using preprocessor directives.

---