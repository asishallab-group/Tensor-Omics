# Makefile -- build the Tensor Omics manual.
#
#   make          build tox_manual.pdf (default)
#   make watch    continuous rebuild + preview (latexmk -pvc)
#   make lint     sanity-check the style package with chktex
#   make clean    remove LaTeX aux files (keeps the PDF)
#   make distclean  remove aux files AND the PDF

MAIN    = tox_manual
LATEXMK = latexmk

.PHONY: all watch lint clean distclean

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex toxmanual.sty $(wildcard recipes/*.tex) $(wildcard shared/*.tex) references.bib
	$(LATEXMK) $(MAIN).tex

watch:
	$(LATEXMK) -pvc $(MAIN).tex

lint:
	chktex toxmanual.sty

clean:
	$(LATEXMK) -c $(MAIN).tex
	@# belt-and-suspenders: mirror misc/clean.sh's aux-extension sweep
	@find . -type f -regextype posix-egrep \
	  -regex '.*(-blx\.bib|\.nav|\.snm|\.aux|\.bbl|\.blg|\.lof|\.log|\.lot|\.run\.xml|\.synctex\.gz|\.toc|\.dvi|\.loa|\.out|\.bcf|\.fdb_latexmk|\.fls|\.ist|\.glo|\.acn|\.vrb)' \
	  -exec rm -f {} +

distclean:
	$(LATEXMK) -C $(MAIN).tex
