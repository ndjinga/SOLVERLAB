#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
example of launch:

>>> matix shell
>>> find . -name "xsal*py"
>>> $PACKAGESPY_SRC_DIR/pythonAppliMatix/salomepy/test/xsalomesessionStudiesTest.py
>>>
>>> AllTestLauncher xsalomesessionStudiesTest.py
>>> import AllTestLauncher as ATL
>>> ATL.runFromEnvVar("PACKAGESPY_ROOT_DIR", "xsalomesessionStudiesTest.py")
"""

import os
import unittest
import salomepy.xsalomesession as XSS
import salomepy.utilsWorkdir as UTW

runningDir = UTW.getWorkdirDefault("TESTS")  

try:
  import salome
  isSalome = True #only user eye: too dangerous for salome desktop
except:
  isSalome = False
  print("WARNING: xsalomesessionStutiesTest: can't import salome: no tests")
  
verbose = False
myDir, myName =  os.path.split(__file__)

class TestCase(unittest.TestCase):
  
  _name = os.path.join(runningDir, "initialActiveStudy.hdf")
  _nameDump = os.path.join(myDir, "simpleBoxAndMeshDump.py")
  _nameHdf = os.path.join(myDir, "simpleBoxAndMesh.hdf")
  _isSalomeAlreadyLaunched = None
  _initialActiveStudy = None
  _activeStudy = None
  _hasDesktopAndAlreadyLaunched = None
  
  def test000(self): #first test
    self.assertTrue(isSalome)
    res = XSS.isSalomeAlreadyLaunched()
    if verbose: print("XSS.isSalomeAlreadyLaunched()",XSS.isSalomeAlreadyLaunched())
    self.assertEqual(res in [True,False], True)
    TestCase._isSalomeAlreadyLaunched = res
    desktop = XSS.getDesktop()
    TestCase._hasDesktopAndAlreadyLaunched = (res and (desktop != None))

  def test010(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    TestCase._initialActiveStudy = XSS.getActiveStudy()
    if verbose: print("XSS.getActiveStudy()",TestCase._initialActiveStudy)

  def test015(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    if os.path.isfile(TestCase._name): os.remove(TestCase._name)
    if TestCase._initialActiveStudy != None:
      # not using multifile save mode
      multifile = False
      if verbose: print("WARNING: delete existing user initial study, save as",TestCase._name)
      salome.myStudyManager.SaveAs(TestCase._name, TestCase._initialActiveStudy, multifile)
      self.assertEqual(os.path.isfile(TestCase._name), True)
    TestCase._activeStudy = XSS.createNewStudy() #only one, delete previous study
    self.assertNotEqual(TestCase._activeStudy, None)

  def test020(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    res = XSS.isSalomeAlreadyLaunched()
    self.assertEqual(res, True)
    self.assertEqual(TestCase._activeStudy, XSS.getActiveStudy())
    self.assertEqual(TestCase._activeStudy, salome.myStudy)

  def test030(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    idstudy1 = id(salome.myStudy)
    ok, why = XSS.openHdf("inexistingFile.tmp")
    self.assertEqual(ok, False)
    ok, why = XSS.openHdf(TestCase._nameHdf)
    self.assertEqual(ok, True)
    self.assertEqual(XSS.getActiveStudy(), salome.myStudy)
    self.assertNotEqual(idstudy1, id(salome.myStudy))
    
  def test110(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    sgeom = salome.myStudy.FindObject("Geometry")
    smesh = salome.myStudy.FindObject("Mesh")
    self.assertEqual(salome.myStudy.GetObjectPath(salome.myStudy.FindObject("Box_1")), "/Geometry/Box_1")
    self.assertEqual(salome.myStudy.GetObjectPath(salome.myStudy.FindObject("Mesh_1")), "/Mesh/Mesh_1")
    objs = []
    objs.append(salome.myStudy.FindObject("Geometry"))
    objs.append(salome.myStudy.FindObject("Mesh"))
    objs.append(salome.myStudy.FindObjectByPath("/Geometry/Box_1"))
    objs.append(salome.myStudy.FindObject("Box_1"))
    objs.append(salome.myStudy.FindObjectByPath("/Mesh/Mesh_1"))
    objs.append(salome.myStudy.FindObject("Mesh_1"))
    #print "objs",objs

    # iterate through objects of the data tree with child iterator
    Verbose = False
    res = XSS.getGeomSobj(Verbose)
    XSS.getOrLoadSobj(res)
    self.assertEqual('/Geometry/Box_1' in XSS.getSobjPath(res), True)
    
    res = XSS.getSmeshSobj(Verbose)
    XSS.getOrLoadSobj(res)
    self.assertEqual('/Mesh/Mesh_1' in XSS.getSobjPath(res), True)
    
    self.assertNotEqual(XSS.getOrLoadSobj(salome.myStudy.FindObject("Box_1")), None)
    self.assertNotEqual(XSS.getOrLoadSobj(salome.myStudy.FindObject("Mesh_1")), None)
    
    """
    matix shell
    ./SOURCES/PACKAGESPY/pythonAppliMatix/salomepy/test/xsalomesessionStudiesTest.py
    $PACKAGESPY_SRC_DIR/pythonAppliMatix/salomepy/test/xsalomesessionStudiesTest.py

    AllTestLauncher xsalomesessionStudiesTest.py
    import AllTestLauncher as ATL
    ATL.runFromEnvVar("PACKAGESPY_ROOT_DIR", "xsalomesessionStudiesTest.py")


    from salome.geom import geomBuilder
    geompy = geomBuilder.New()
    from salome.smesh import smeshBuilder
    smeshpy = smeshBuilder.New()
    
    salome.salome_init(salome.myStudy)
    g = salome.lcc.FindOrLoadComponent("FactoryServer", "GEOM")
    s = salome.lcc.FindOrLoadComponent("FactoryServer", "SMESH")
    print("Component",g,s)
    g.init_geom(salome.myStudy)
    XSS.getGeomGUI()
    """

  def xtest120(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    desktop = XSS.getDesktop()
    if desktop == None:
      self.assertEqual(XSS.getSelectedCount(), None)
    else:
      self.assertEqual(XSS.getSelectedCount()>=0, True)
    
  def xtest130(self):
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    desktop = XSS.getDesktop()
    if desktop == None:
      self.assertEqual(XSS.getAllSelected(), None)
    else:
      self.assertEqual(len(XSS.getAllSelected())>=0, True)
    
  def xtest130(self):
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
    
  def xtest999(self): #last test
    #avoid th. SALOME_NamingService.cxx [1379] : Destroy_Directory(): CosNaming::NamingContext::NoEmpty /Study/ is not empty
    #avoid kill previous user study...
    if not isSalome: return
    self.assertFalse(TestCase._hasDesktopAndAlreadyLaunched)
    XSS.closeActiveStudy()
    self.assertEqual(XSS.getActiveStudy(), None)
    for i in dir(salome.myStudyManager): print(i)
    if os.path.isfile(TestCase._name): 
      study = salome.myStudyManager.Open(TestCase._name)
      self.assertEqual(XSS.getActiveStudy(), salome.myStudy)
      self.assertEqual(XSS.getActiveStudy(), study)

if __name__ == '__main__':
  verbose = True   #verbose if human eye
  unittest.main()

