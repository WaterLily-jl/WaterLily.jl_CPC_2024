MAIN=main
OUT_DIR=build

.SUFFIXES:
.SUFFIXES: .bib .pdf .tex
.PHONY: clean

run: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).bbl $(MAIN).tex
	pdflatex --output-directory=$(OUT_DIR) $(MAIN).tex -draftmode
	pdflatex --output-directory=$(OUT_DIR) $(MAIN).tex
	cp $(OUT_DIR)/$(MAIN).pdf .

$(MAIN).bbl: $(MAIN).aux
	bibtex $(OUT_DIR)/$(MAIN)

$(MAIN).aux: $(MAIN).bib
	pdflatex --output-directory=$(OUT_DIR) $(MAIN).tex -draftmode
	pdflatex --output-directory=$(OUT_DIR) $(MAIN).tex -draftmode

clean:
	rm -rf *.aux *.lof *.log *.lot *.toc *.bbl *.blg *.pdf *.out build

 $(info $(shell mkdir -p $(OUT_DIR)))