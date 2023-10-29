
### SOLVERLABGUI IN MATIX

Compilation installation in SALOME/MATIX environment with `CMAKE` usage.

See CMakeLists.txt etc.


#### git create branch 

```
git checkout -b cvw_M26_S850_SALOME

new addings
  README_MATIX.md

addings for SALOME modules (as DART for example)
  bin
  resources
  src
  CMakeLists.txt
```

#### compilation example with sat

```
ADIR=/volatile2/wambeke/TULEAP_MATIX
TRG=RELEASED_MATIX_V26_850_COS73
cd ${ADIR}
alias ssat=${ADIR}/SAT/sat

rm -rf MATIX_26-CO7/BUILD/SOLVERLABGUI MATIX_26-CO7/INSTALL/SOLVERLABGUI

ssat compile $TRG -p SOLVERLABGUI --clean_install
ssat compile $TRG -p MATIX_PROFILE --clean_install
cp -r MATIX_26-CO7/SOURCES/SOLVERLABGUI/*py MATIX_26-CO7/INSTALL/SOLVERLABGUI

ssat launcher $TRG
ssat launcher $TRG --use_mesa -n matix_mesa

${ADIR}/MATIX_26-CO7/matix -k
```