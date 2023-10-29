#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
example of launch:

>>> matix shell
>>> find . -name "xsal*py"
>>> $PACKAGESPY_SRC_DIR/pythonAppliMatix/salomepy/test/xsalomesessionTest.py
>>>
>>> import AllTestLauncher as ATL
>>> ATL.runFromEnvVar("PACKAGESPY_ROOT_DIR", "xsalomesessionTest.py")
"""

import os
import unittest

verbose = False
import salomepy.xsalomesession as XSS

try:
  import salome
  isSalome = True
except:
  isSalome = False
  print("WARNING: xsalomesessionTest: can't import salome: no tests")

class TestCase(unittest.TestCase):
  
  _isSalomeAlreadyLaunched = None
  _initialActiveStudy = None
  _hasDesktopAndAlreadyLaunched = None
  
  def test_000(self):
    self.assertTrue(isSalome)
    res = XSS.isSalomeAlreadyLaunched()
    if verbose: print("XSS.isSalomeAlreadyLaunched()",XSS.isSalomeAlreadyLaunched())
    self.assertEqual(res in [True,False], True)
    TestCase._isSalomeAlreadyLaunched = res
    desktop = XSS.getDesktop()
    TestCase._hasDesktopAndAlreadyLaunched = (res and (desktop != None))

  def test_010(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    TestCase._initialActiveStudy = XSS.getActiveStudy()
    if verbose: print("XSS.getActiveStudy()",TestCase._initialActiveStudy)
    
  def test_020(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    if TestCase._isSalomeAlreadyLaunched == False:
      XSS.getXSalomeSession() #create a default study...
    res = XSS.isSalomeAlreadyLaunched()
    self.assertEqual(res, True)

  def test_030(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    if TestCase._initialActiveStudy != None:
      print("WARNING: xsalomesessionTest: do not test to create study deleting user initial study")
      return
    newStudy = XSS.createNewStudy()
    if verbose: print("newStudy", newStudy)
    self.assertNotEqual(newStudy, None)
    
  def test_110(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    from PyQt5.QtGui import QMainWindow
    desktop = XSS.getDesktop()
    self.assertEqual(type(desktop) in [QMainWindow, type(None)], True)
    
  def test_120(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    desktop = XSS.getDesktop()
    if desktop == None:
      self.assertEqual(XSS.getSelectedCount(), None)
    else:
      self.assertEqual(XSS.getSelectedCount()>=0, True)
    
  def test130(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    desktop = XSS.getDesktop()
    if desktop == None:
      self.assertEqual(XSS.getAllSelected(), None)
    else:
      self.assertEqual(len(XSS.getAllSelected())>=0, True)
    
  def test130(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    desktop = XSS.getDesktop()
    if desktop == None:
      self.assertEqual(XSS.getSelected(), None)
    else:
      res=XSS.getSelected(0)
      if XSS.getSelectedCount() == 0:
         self.assertEqual(res, None)
      else:
         self.assertEqual(type(res), str) #entry as '0:1:2' for example
    
  def test_999(self):
    #avoid th. SALOME_NamingService.cxx [1379] : Destroy_Directory(): CosNaming::NamingContext::NoEmpty /Study/ is not empty
    #avoid kill previous user study...
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    if TestCase._initialActiveStudy != None:
      print("WARNING: xsalomesessionTest: do not test to close an user initial study")
      return
    XSS.closeActiveStudy()
    self.assertEqual(XSS.getActiveStudy(), None)


if __name__ == '__main__':
  verbose = True   #verbose if human eye
  unittest.main()

