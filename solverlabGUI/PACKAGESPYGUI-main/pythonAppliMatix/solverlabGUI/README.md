

### User documentation

The solverlab GUI is python3 package which currently launch solverlab CODE.
Have to set a python3 environment.

Solverlab CODE is compiled C code.

Documentation solverlab GUI, see `.../doc/build/html/index.html`

Documentation solverlab CODE, see `.../doc/src/solverlabDocuments/20140804_solverlab_manual.pdf`


### Release information

see .../sandbox/TODO_*.txt


### Developer information

To create archives linux `.zip` `.tgz` or windows `.7z`.
see `.../sandbox/README_create_solverlab_zip_tgz_7z.md`

To create python self extracting package by pymakeself (for future not done yet 2020).
See `.../solverlabGUI/sandbox/pymakeself/README_pymakeself.txt`

Math tex expressions information.
see `https://matplotlib.org/users/mathtext.html#mathtext-tutorial`



### Set python environment

To get python3 PyQt5 correct environment, use conda/miniconda.
see `https://docs.conda.io/en/latest/miniconda.html`.
see `.../solverlabGUI/sandbox/conda/README_pySolverlabGUI3.md`

This is a conda/miniconda example.

```
export CONDADIR="/volatile/common/miniconda3"
export PATH=$CONDADIR/bin:$PATH
conda activate pySolverlabGUI3.7   # python3 pyqt5 numpy pandas matplotlib...
```


### Launch solverlab GUI


Launch main script is `.../solverlabGUI/solverlabGUI`

```
cd .../solverlabGUI

./solverlabGUI --help   # get help
./solverlabGUI --doc    # get documentation

./solverlabGUI --gui    # launch GUI
```

### Launch solverlab GUI simple Tests


```
cd .../solverlabGUI
./AllTestLauncher.sh
```

### Launch solverlab CODE integration tests


```
cd .../solverlabGUI/solverlabCodes
cat README.txt
./AllTestLauncher.sh
```