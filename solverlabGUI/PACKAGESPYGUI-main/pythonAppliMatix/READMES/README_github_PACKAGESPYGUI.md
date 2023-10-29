

## Github PACKAGESPYGUI

Github [PACKAGESPYGUI](https://github.com/ndjinga/PACKAGESPYGUI)
is a **fork** (halas) from CEA tuleap SALOME/MATIX
[PACKAGESPY](https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy).

**warning:**
this is tricky way forking Tuleap PACKAGESPY to Github PACKAGESPYGUI as master
because they are *not* common history shared git repositories.


### Fork Tuleap to github

```
# all this tricky way is in temporary directory
mkdir tmp
cd tmp
# mandatory name PACKAGESPY_TULEAP as precaution
rm -rf PACKAGESPY_TULEAP
git clone https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy PACKAGESPY_TULEAP
# set your branch for PACKAGESPY_TULEAP
cd PACKAGESPY_TULEAP
git checkout cv_M30_2022

./configure_solverlabGUI.sh --help
# choose one of these
#./configure_solverlabGUI.sh solverlab_standalone
#./configure_solverlabGUI.sh solverlab_in_matix
./configure_solverlabGUI.sh solverlab_in_salome

# push new fork of PACKAGESPY in https://github.com/ndjinga/PACKAGESPYGUI as master
# tricky way because that is not common shared git repository
cd ..
git clone https://github.com/ndjinga/PACKAGESPYGUI
rm -rf PACKAGESPYGUI/*
cp -r PACKAGESPY_TULEAP/* PACKAGESPYGUI/.
cd PACKAGESPYGUI
git st
git add -A
git commit -m "fork from $(tail -n 1 git_history.tmp)"
git push
```
