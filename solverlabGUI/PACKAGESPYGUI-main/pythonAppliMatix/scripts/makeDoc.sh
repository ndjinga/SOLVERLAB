#!/bin/bash

echo "fromFilePath="`pwd`/$0
cd doc
make html
make latexpdf

echo "lire la doc html: firefox ./doc/_build/html/index.html"
echo "lire la doc pdf : evince  ./doc/_build/latex/packagespy.pdf"

