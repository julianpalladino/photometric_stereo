.PHONY: clean open

PDF=informe.pdf

$(PDF): informe.tex $(wildcard *.tex)
	latexmk -pdf $<

clean:
	latexmk -C

open:
	xdg-open $(PDF)
