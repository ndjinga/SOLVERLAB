on wambeke@is246206

f_init_conda
  # conda environments:
  #
  base                  *  /export/home/catA/miniconda
  py3qt5                   /export/home/catA/miniconda/envs/py3qt5


@is246206.../PACKAGESPY(cv_M30_2022)>conda activate py3qt5
@is246206.../PACKAGESPY(cv_M30_2022)>which python
/export/home/catA/miniconda/envs/py3qt5/bin/python
@is246206.../PACKAGESPY(cv_M30_2022)>cd /export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC

echo $PYTHONPATH

export PYTHONPATH=/export/home/catA/wambeke/SALOME/SALOME-9.9.0-native-FD34-SRC/INSTALL/PACKAGESPY/pythonAppliMatix

./INSTALL/PACKAGESPY/pythonAppliMatix/amitexpy/guiAmt/test/test_510_treeXmlAmt.py

.
----------------------------------------------------------------------
Ran 2 tests in 9.707s

OK



export aDir=/home/christian/Documents/tmp
#or
export aDir=/export/home/matix/tmp_packagespy/packagespy
cd ~

#export PYTHONPATH=$aDir/pythonAppliMatix
#echo $PYTHONPATH

bash
cd ~
export PYTHONPATH=$aDir/pythonAppliMatix
echo $PYTHONPATH
#erreur!!
python
import pythonAppliMatix as TMP
for i in dir(TMP):
  print i

exit()
exit


bash
cd ~
export PYTHONPATH=$aDir/pythonAppliMatix
echo $PYTHONPATH
python
import xyzpy as TMP
for i in dir(TMP):
  print i

import xyzpy.test as TMP

exit()
exit


bash
cd ~
export PYTHONPATH=$aDir/pythonAppliMatix
echo $PYTHONPATH
python
import xyzpy.modelXyz as TMP
for i in dir(TMP):
  print i

exit()
exit


bash
cd ~
export PYTHONPATH=$aDir/pythonAppliMatix
echo $PYTHONPATH
python
import xyzpy.test as TMP
for i in dir(TMP):
  print i

from xyzpy.test import *
exit()
exit


bash
/export/home/wambeke/2014/MATIXP_V2_740_CV_CO6.4_64/matix_S740 shell
export PYTHONPATH=$aDir/pythonAppliMatix
echo $PYTHONPATH
cd $aDir/pythonAppliMatix/xyzpy/test
baseXyzTest.py



python
import xyzpy.AllTestLauncher as TT
TT.run()

import xyzpy.test.baseXyzTest as TT
TT.run() #do not work

exit()



#debug
ssh -Y matix@is206786
export aDir=/export/home/matix/tmp_packagespy/packagespy
bash
/export/home/wambeke/2014/MATIXP_V2_740_CV_CO6.4_64/matix_S740 shell
cd /export/home/matix
export PYTHONPATH=$aDir/pythonAppliMatix:$PYTHONPATH
echo $PYTHONPATH
tree $PYTHONPATH




python
import xyzpy.test.baseXyzTest


import xyzpy.AllTestLauncher
xyzpy.AllTestLauncher.run()

exit()
exit

cd $aDir/pythonAppliMatix/xyzpy/guiXyz/test
