#!/bin/bash
find . -type f -regextype posix-egrep \
  -regex '.*(-blx\.bib|\.nav|\.snm|\.aux|\.bbl|\.blg|\.lof|\.log|\.lot|\.run\.xml|\.synctex\.gz|\.toc|\.dvi|\.loa|\.out|\.bcf|\.fdb_latexmk|\.fls|\.ist|\.glo|\.acn|\.vrb)' \
  -exec rm -f \{\} +;
