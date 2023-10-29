#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
singleton on QtWidgets.QApplication
"""

import sys
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QT_VERSION_STR
# cvw TODO import xyzpy.loggingXyz as LOG

__app__ = [None]
verbose = False

# logger = LOG.getLogger()

# avoid in 'make html' typeError: iter() returned non-iterator of type 'NoneType'
try:
  versionQt = [int(i) for i in QT_VERSION_STR.split('.')]
except:
  versionQt = [5, 0, 0]

versionPython = sys.version_info

if versionPython > (3, 0,) and versionPython < (3, 6,):
  if versionQt >= (5, 0,) and versionQt < (5, 9,):
    print("""\
WARNING : QApplication (may be) have segmentation fault
see https://github.com/pyinstaller/pyinstaller/issues/2154
   seg fault on app = QApplication(sys.argv)
   the issue is gone with PyInstaller-3.2.1 PyQt5.7.1 Python 3.5.2
   (salome 9 is PyQt5.9 Python 3.6.0, ouf)
""" )

def OnceQApplication(*args):
  # print("WARNING : !!! only one QApplication have to be done")
  if __app__[0] is None:
    if args == ():
      __app__[0] = QApplication([])
    else:
      __app__[0] = QApplication(*args)
    try:
      app__[0].setStyle("Fusion") #['Oxygen', 'Windows', 'Motif', 'CDE', 'Plastique', 'GTK+', 'Cleanlooks']
    except: # may be not existing
      pass
    if verbose: print("WARNING : !!! QApplication done NOW")
  else:
    if args != ():
      pass
      # cvw TODO logger.warning("QApplication done yet: args are not considered: %s" % args)
      # print("WARNING : !!! QApplication done YET but args are not considered: %s" % args)
    else:
      pass
    if verbose: print("WARNING : !!! QApplication done YET")
  return __app__[0]


if __name__ == '__main__':
  import sys
  app1=OnceQApplication(sys.argv)
  print(id(app1))
  app2=OnceQApplication(sys.argv)
  print(id(app2))
  app3=OnceQApplication()
  print(id(app3))
