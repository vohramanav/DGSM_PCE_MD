#!/bin/sh

pdflatex outline.tex
pdflatex outline.tex
rm *.aux *.log
