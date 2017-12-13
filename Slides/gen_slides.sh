#!/bin/sh

pdflatex dgsm_pce.tex
pdflatex dgsm_pce.tex
rm *.aux *.log *.out *.toc *.nav *.snm
