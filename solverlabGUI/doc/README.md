### MAKE documentation (LINUX)

- Needs python sphinx-build to make doc html
- solverlabGUI.pdf and html's are pushed in git
  to get up-to-date documentation directly from git clone
- Compilation sphinx have to be *done before* MATIX-sat compilation.

```
# may be to get python/sphinx prerequisites
matix context
# or
conda activate py3

cd solverlabGUI

# to make doc html apidoc on commands dir ok
export PYTHONPATH=$(pwd)

cd doc
rm -rf buid/*
make html
firefox build/html/index.html &

# to make doc pdf
# needs texlive up to date (done on machines lgls for x86_64-linux)
# https://www.tug.org/texlive/quickinstall.html

# an example of set PATH's
export TEX_ROOT_PATH=/data/tmplgls/wambeke/share/texlive/2017
export INFOPATH=${TEX_ROOT_PATH}/texmf-dist/doc/info
export MANPATH=${TEX_ROOT_PATH}/texmf-dist/doc/man
export PATH=${TEX_ROOT_PATH}/bin/x86_64-linux:${PATH}

cd doc
make latexpdf
# makefile copy build/latex/solverlabGUI.pdf in src/solverlabDocuments
evince build/latex/solverlabGUI.pdf &
evince src/solverlabDocuments/solverlabGUI.pdf &
```
