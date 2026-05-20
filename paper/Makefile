LATEX = pdflatex
BIBTEX = bibtex

.PHONY: all check clean

all: main.pdf

main.pdf: main.tex refs.bib llncs.cls splncs04.bst
	$(LATEX) main
	$(BIBTEX) main
	$(LATEX) main
	$(LATEX) main

check: main.pdf
	pdfinfo main.pdf
	pdffonts main.pdf

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.pdf main.toc
