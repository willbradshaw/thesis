#!/bin/bash
# A script to compile the PhD Thesis - Krishna Kumar 
# Distributed under GPLv2.0 License

TEXDIR="/usr/local/texlive/2018/bin/x86_64-linux"
TEXPATH="${TEXDIR}/pdflatex"
BIBPATH="${TEXDIR}/biber"
GLOPATH="${TEXDIR}/makeglossaries"
INDPATH="${TEXDIR}/makeindex"
compile="compile";
clean="clean";

if test -z "$2"
then
if [ $1 = $clean ]; then
	echo "Cleaning please wait ..."
	rm -f *~
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.d
	rm -rf *.fls
	rm -rf *.ilg
	rm -rf *.ind
	rm -rf *.toc*
	rm -rf *.lot*
	rm -rf *.lof*
	rm -rf *.log
	rm -rf *.idx
	rm -rf *.out*
	rm -rf *.nlo
	rm -rf *.nls
	rm -rf *.bcf
	rm -rf *.run.xml
	rm -rf *.tex.bak
	rm -rf *.glg
	rm -rf *.glo
	rm -rf *.gls
	rm -rf *.glsdefs
	rm -rf *.ist
	rm -rf $filename.pdf
	rm -rf $filename.ps
	rm -rf $filename.dvi
	rm -rf *#* 
	echo "Cleaning complete!"
	exit
else
	echo "Shell script for compiling the PhD Thesis"
	echo "Usage: sh ./compile-thesis.sh [OPTIONS] [filename]"
	echo "[option]  compile: Compiles the PhD Thesis"
	echo "[option]  clean: removes temporary files no filename required"
	exit
fi
fi

filename=$2;

if [ $1 = $clean ]; then
	echo "Cleaning please wait ..."
	rm -f *~
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.d
	rm -rf *.fls
	rm -rf *.ilg
	rm -rf *.ind
	rm -rf *.toc*
	rm -rf *.lot*
	rm -rf *.lof*
	rm -rf *.log
	rm -rf *.idx
	rm -rf *.out*
	rm -rf *.nlo
	rm -rf *.nls
	rm -rf *.bcf
	rm -rf *.glg
	rm -rf *.glo
	rm -rf *.gls
	rm -rf *.glsdefs
	rm -rf *.ist
	rm -rf *.run.xml
	rm -rf *.tex.bak
	rm -rf $filename.pdf
	rm -rf $filename.ps
	rm -rf $filename.dvi
	rm -rf *#* 
	echo "Cleaning complete!"
	exit
elif [ $1 = $compile ]; then
	echo "Compiling your PhD Thesis...please wait...!"
	$TEXPATH -interaction=nonstopmode $filename.tex
	$BIBPATH $filename
	$GLOPATH $filename
	$INDPATH $filename.aux
	$INDPATH $filename.idx
	$INDPATH $filename.nlo -s nomencl.ist -o $filename.nls
	$TEXPATH -interaction=nonstopmode $filename.tex
	$INDPATH $filename.nlo -s nomencl.ist -o $filename.nls
	$GLOPATH $filename
	$TEXPATH -interaction=nonstopmode $filename.tex
    echo "Cleaning auxiliary files..."
	rm -rf */*.aux
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.d
	rm -rf *.fls
	rm -rf *.ilg
	rm -rf *.ind
	rm -rf *.toc*
	rm -rf *.lot*
	rm -rf *.lof*
	rm -rf *.log
	rm -rf *.idx
	rm -rf *.out*
	rm -rf *.nlo
	rm -rf *.nls
	rm -rf *.bcf
	rm -rf *.glg
	rm -rf *.glo
	rm -rf *.gls
	rm -rf *.glsdefs
	rm -rf *.ist
	rm -rf *.run.xml
	rm -rf *.tex.bak
	echo "Success!"
	exit
fi


if test -z "$3"
then
	exit
fi
