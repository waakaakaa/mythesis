FileName=论文模版示例
PDFPara=-p a4
TEXPara=--src-specials --synctex=-1

$(FileName).pdf : $(FileName).dvi
	dvipdfm $(PDFPara) $(FileName).dvi
	-del *.dvi	

$(FileName).dvi : $(FileName).tex
	latex $(TEXPara) $(FileName)
	makeindex $(FileName).idx
	bibtex $(FileName) 
	latex $(TEXPara) $(FileName)
	latex $(TEXPara) $(FileName)

clear :
	-del $(FileName).dvi
	-del $(FileName).aux
	-del $(FileName).log
	-del $(FileName).pdf
	-del $(FileName).toc
	-del $(FileName).idx
	-del $(FileName).ind
	-del $(FileName).out
	-del $(FileName).bbl
	-del $(FileName).lof
	-del $(FileName).lot
	-del $(FileName).ilg
	-del $(FileName).blg
	-del $(FileName).synctex
#	-del bu*.*	

