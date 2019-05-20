########################################################
all:
	pdflatex -interaction nonstopmode qMRI-Example.tex
	-bibtex qMRI-Example.aux
	pdflatex -interaction nonstopmode qMRI-Example.tex
	pdflatex -interaction nonstopmode qMRI-Example.tex

clean:
	rm *.aux *.bbl *.blg *.log *.tex *.bib *.xml
