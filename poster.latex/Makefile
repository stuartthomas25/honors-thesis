.PHONY: all

all : pdf clean

pdf: main.tex $(wildcard imgs/*)
	#latexmk -pdflatex="pdflatex -synctex=1 %O %S" -pdf main.tex
	latexmk -pdf main.tex

.PHONY: clean

clean:
	latexmk -c

