# User-configurable stuff -------------------------------------------
BASE=ramses_ug
LOG=userfiles/tube1d.log
NML=userfiles/tube1d.nml
AWK=gawk
AUTODIR=autolog
# -------------------------------------------------------------------

TEX=$(BASE).tex
ALLTEX=*.tex
PDF=$(BASE).pdf
IND=$(BASE).ind
IDX=$(BASE).idx

pdf: $(ALLTEX) $(IND) autolog.target
	pdflatex $(TEX)

$(IND): $(IDX)
	makeindex $(IDX)

$(IDX): $(ALLTEX) autolog.target
	pdflatex $(TEX)

autolog.target:
	mkdir -p $(AUTODIR)
	scripts/autolog.sh $(LOG) $(NML) $(AUTODIR)

clean:
	rm -f *.aux *.idx *.ilg *.ind *.log *.out *.toc \
		$(AUTODIR)/*.tex $(AUTODIR)/*.log $(AUTODIR)/*.aux

allclean: clean
	rm -f *.pdf
