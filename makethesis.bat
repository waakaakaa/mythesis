latex --src-specials --synctex=-1 main
makeindex main.idx
bibtex main
latex --src-specials --synctex=-1 main
latex --src-specials --synctex=-1 main
dvipdfmx -p a4 main
del *.dvi
del *.aux
del *.log
del *.toc
del *.idx
del *.ind
del *.out
del *.bbl
del *.lof
del *.lot
del *.ilg
del *.blg
del *.synctex