#!/bin/sh

pdflatex DGSM_PCE_MD_paper.tex
pdflatex DGSM_PCE_MD_paper.tex
rm *.aux *.log *.toc
