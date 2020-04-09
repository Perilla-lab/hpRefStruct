PROJNAME=main
VERSION=0.20
GIT_BIBROOT=https://github.com/jrperilla/bibfiles.git
BIB_GITNAME=bibfiles
GIT=git
DETEX=detex

default: pdf

pdf:
	-pdflatex $(PROJNAME)
	-bibtex $(PROJNAME)
	pdflatex $(PROJNAME)
	pdflatex $(PROJNAME)

clean:
	-/bin/rm $(PROJNAME).aux $(PROJNAME).bbl $(PROJNAME).blg $(PROJNAME).log $(PROJNAME).pdf


bibfiles:
	$(GIT) clone $(GIT_BIBROOT)

checkin: 
	$(GIT) push

checkout:
	$(GIT) pull

spell:
	cat *.tex | detex | spell | more

wordcount:
	-echo `$(DETEX) $(PROJNAME).tex | wc -w | tr -d [:space:]` Words

