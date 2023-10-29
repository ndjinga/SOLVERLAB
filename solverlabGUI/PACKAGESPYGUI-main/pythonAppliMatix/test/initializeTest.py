#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
initialize PATH etc... for iradinaGui and/or packagespy unittest test files
"""

"""
https://docs.python.org/2/library/os.html
os.environ mapping is captured the first time the os module is imported, 
typically during Python startup as part of processing site.py. 
Changes to the environment made after this time are not reflected 
in os.environ, except for changes made by modifying os.environ directly.

On some platforms, including FreeBSD and Mac OS X, 
setting environ may cause memory leaks.
"""

import os
import sys
import pprint as PP

def addSysPath(aDir):
  if aDir not in sys.path:
    # get path to packagespy or iradinaGui sources FIRST as prepend
    # Make the src & commands package accessible from all unittest code
    sys.path.insert(0, aDir)
    sys.stderr.write("""\
WARNING : sys.path not set for unittest, fixed for you:
          sys.path prepend '%s'
          sys.path:\n%s\n""" % (aDir, PP.pformat(sys.path)))
    # os.environ PATH is not set...
    sys.stderr.write("INFO    : to fix this message type:\n  'export PATH=%s:${PATH}'\n" % aDir)

# get path to packagespy or iradinaGui sources directory parent
# aDir = os.path.realpath(os.getcwd()) 
# aDir is something as ".../PACKAGESPY/pythonAppliMatix/test" or ".../IRADINAGUI/test"
aDir = os.path.dirname(os.path.realpath(__file__))

if "PACKAGESPY" in aDir.upper():
  # in salome INSTALL packagespy
  addSysPath(os.path.join(aDir.split("pythonAppliMatix")[0], "pythonAppliMatix"))

elif "IRADINAGUI" in aDir.upper():
  # in standalone iradinagui
  addSysPath(aDir.split("test")[0])


