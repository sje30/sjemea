all: sjemea-intro.pdf

sjemea-intro.pdf: sjemea-intro.Rnw sjemea.bib
	Rscript -e 'require(knitr); knit2pdf("sjemea-intro.Rnw")'

clean:
	rm -f waverepo_paper.tex *.blg *.bbl waverepo_paper.pdf


cleanall:
	rm -fr cache sjemea-intro.pdf
