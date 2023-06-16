pdflatex libcKrylov_doc.tex
bibtex libcKrylov_doc.aux
pdflatex libcKrylov_doc.tex
pdflatex libcKrylov_doc.tex


rm -rf libcKrylov_doc.aux
rm -rf libcKrylov_doc.dvi
rm -rf libcKrylov_doc.log
rm -rf libcKrylov_doc.toc
rm -rf libcKrylov_doc.bbl
rm -rf libcKrylov_doc.blg
rm -rf libcKrylov_doc.out
rm *~

evince libcKrylov_doc.pdf &
