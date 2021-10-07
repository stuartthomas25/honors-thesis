.PHONY: all

all : pdf clean

pdf: main.tex $(wildcard imgs/*)
	latexmk -pdf main.tex

.PHONY: clean

clean:
	latexmk -c main.tex

