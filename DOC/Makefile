
SHELL = /bin/sh

TEXHOME = /salsify/packages/teTeX/bin/
TEXHOME = /usr/bin/
TEXHOME = /usr/cluster/contrib/texlive2018/bin/x86_64-linux/pdf
TEXHOME = /cm/shared/apps/texLive/2025/bin/x86_64-linux/pdf

all:kcarta1.22.pdf kcarta1.12.pdf ReadMeFirst.pdf

kcarta1.22.pdf: kcarta1.22.tex ASL_defs.tex
	$(TEXHOME)latex kcarta1.22.tex && \
	$(TEXHOME)latex kcarta1.22.tex && \
	bibtex kcarta1.22 && \
	$(TEXHOME)latex kcarta1.22.tex

kcarta1.12.pdf: kcarta1.12.tex ASL_defs.tex
	$(TEXHOME)latex kcarta1.12.tex && \
	$(TEXHOME)latex kcarta1.12.tex 

ReadMeFirst.pdf: ReadMeFirst.tex ASL_defs.tex
	$(TEXHOME)latex ReadMeFirst.tex && \
	$(TEXHOME)latex ReadMeFirst.tex

clean:
	rm *.aux *.lof *.dvi *.log *.lot *.toc *.ps *.pdf *~ 2> /dev/null || true


