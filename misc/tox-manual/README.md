# Tensor Omics manual — `toxmanual.sty`

This directory holds the **new, recipe-style TOX user manual** (issue
[#147](https://github.com/asishallab-group/Tensor-Omics/issues/147)) and the
LaTeX style package that makes writing it easy and consistent.

The manual is a **cookbook of analysis "recipes"**, styled after Bioconductor
vignettes / the edgeR User's Guide. It documents reproducible *workflows*, not
API or source code. It is deliberately separate from the older API-reference
manual on the `general_documentation` branch.

## Layout

```
misc/tox-manual/
  toxmanual.sty        the style package (all commands/environments)
  tox_manual.tex       master document — \input each recipe here
  recipes/
    _template.tex      copy this to start a new recipe
    outlier_detection.tex  divergent-gene analysis
    trajectory_contributions.tex  signed contribution analysis
    snd_detection.tex      SND analysis
  shared/
    tox_dataset.tex             creating a TOX dataset
    normalization.tex           the normalization sub-recipe
    family_centroids.tex        the family-centroid sub-recipe
    gene_family_detection.tex   families and orthologs from OrthoFinder
  references.bib       starter bibliography
  latexmkrc  Makefile  build tooling (pinned to pdfLaTeX)
```

## Building

Requires a full TeX Live / MiKTeX install (provides `tcolorbox`, `newtx`,
`scrreprt`, `fontawesome5`, `cleveref`, `latexmk`).

```sh
make           # build tox_manual.pdf
make watch     # rebuild + preview on every save
make lint      # chktex sanity-check of toxmanual.sty
make clean     # remove aux files (keep the PDF)
```

## Writing a recipe

1. Copy `recipes/_template.tex` to `recipes/<slug>.tex` (or
   `shared/<slug>.tex` for a reusable sub-recipe).
2. Give it a unique label: `[label=rec:<slug>]` for analyses,
   `[label=sub:<slug>]` for shared sub-recipes.
3. Fill in every field. Prefer writing "None" over omitting a field.
4. `\input` your file from `tox_manual.tex` under the right `\part`.

Keep it **simple, safe, clear, concise**. Link to full scripts with
`\examplescript`; show only the fragments a user must inspect or edit.

## Command / environment reference

### Document skeleton
| Command | Purpose |
| --- | --- |
| `\toxtitlepage` | Branded title page. Override with `\toxtitle{}`, `\toxsubtitle{}`, `\toxauthors{}`, `\toxinstitution{}`. |
| `\begin{recipe}[label=rec:x]{Title} … \end{recipe}` | One recipe = one chapter. `label` enables cross-references. |

### Recipe card (top of each recipe)
| Command / environment | Purpose |
| --- | --- |
| `\recipepurpose{…}` | One–two-sentence purpose box. |
| `\nomaintainer` | Marks the recipe as team-owned. **Use this by default.** Prints the name set by `\toxteammaintainer{}` (default: Tensor Omics Developers). |
| `\recipemaintainer{Name}` | Names an individual owner instead. Only when one person really owns the workflow. |
| `\begin{requiredinputs} \item … \end{requiredinputs}` | Data the user must have. |
| `\begin{prerequisites} \dependson{sub:x} \end{prerequisites}` | Depends-on box; each `\dependson` links to another recipe by its label and prints its title + chapter number. |

### Body sections
| Environment | Purpose |
| --- | --- |
| `\begin{workflow} \step{…} \end{workflow}` | Numbered step-by-step procedure. |
| `\begin{parameters} \param{name}{controls}{when/which} \end{parameters}` | Colored parameter table. |
| `\begin{expectedresults} … \end{expectedresults}` | Intermediate results to expect. |
| `\begin{pitfall}[optional title] … \end{pitfall}` | Common pitfall (repeatable). |
| `\begin{interpretation} … \end{interpretation}` | How to read the final output. |
| `\begin{references} \reference{…} \end{references}` | Publications (or use `\cite` + `references.bib`). |

### Code & links
| Command / environment | Purpose |
| --- | --- |
| `\examplescript{path/in/repo}{caption}` | Boxed link to a full script in the repo. |
| `\begin{codefragment}[language=…,caption=…] … \end{codefragment}` | Styled listing for small editable fragments. |
| `\toxrepourl{…}` | Change the repo base URL used by `\examplescript` (default: GitHub `blob/main/`). |
| `\file{}` `\cmd{}` `\code{}` `\optn{}` | Inline monospace helpers. |

### Generic admonitions (usable anywhere)
`\begin{note}…\end{note}` · `\begin{tip}…\end{tip}` · `\begin{warning}…\end{warning}`

### Cross-references
`\seerecipe{rec:x}` prints the recipe's title and chapter number inline;
`\dependson{rec:x}` does the same as a bullet inside `prerequisites`.

## Branding

All colors live at the top of `toxmanual.sty` (`toxblue` = the repo's link blue,
`RGB 0,0,238`). Change them in one place to restyle the whole manual. Body font
is `newtx` (Times-like); there is no logo asset yet.
