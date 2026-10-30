all: report.pdf

report.pdf: report.tex rates_table.tex scales_table.tex ratios_table.tex mass_ratio_table.tex $(wildcard figures/*.pdf)
	pdflatex -interaction=nonstopmode -halt-on-error report.tex
	pdflatex -interaction=nonstopmode -halt-on-error report.tex
