pdflatex libcKrylov.tex
bibtex libcKrylov.aux
pdflatex libcKrylov.tex
pdflatex libcKrylov.tex
#dvipdfm libcKrylov.dvi

rm -rf *.aux
rm -rf *.dvi
rm -rf *.log
rm -rf *.toc
rm -rf *.bbl
rm -rf *.blg
rm -rf *.out
rm -rf *.spl
rm -rf *.gz
rm -rf *.fff
rm -rf *.lof
rm -rf *~

evince libcKrylov.pdf &
 
