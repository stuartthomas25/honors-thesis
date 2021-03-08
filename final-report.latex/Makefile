.PHONY: all

all : pdf clean

pdf: main.tex $(wildcard imgs/*)
	latexmk -pdf

.PHONY: clean

clean:
	latexmk -c

