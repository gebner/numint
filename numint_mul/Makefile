.PHONY: all clean

all: MehrdimensionaleQuadratur.pdf

MehrdimensionaleQuadratur.pdf: MehrdimensionaleQuadratur.tex

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.eps *.pdf *.log *.aux
