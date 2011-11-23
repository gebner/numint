.PHONY: all clean

all: notes.pdf

notes.pdf: notes.tex

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.eps *.pdf *.log *.aux
