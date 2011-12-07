.PHONY: all clean

all: NumerischeQuadratur.pdf

NumerischeQuadratur.pdf: NumerischeQuadratur.tex

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.eps *.pdf *.log *.aux
