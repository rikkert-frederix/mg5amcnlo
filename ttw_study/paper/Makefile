# Numerical tables and figures are generated separately from immutable reports.
.PHONY: all
all: main.pdf

main.pdf: $(wildcard *.tex) jheppub.sty $(wildcard ../figures/*.pdf) $(wildcard ../figures/main_onshell_v3_coarse_v1/*.pdf)
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
