#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""\
tests for dataFromFileXyz
"""

import os
import unittest
import pprint as PP

try:
  import pandas as pd
  pandasOk = True
except:
  print("ERROR: no pandas library, no tests")
  pandasOk = False


if pandasOk:

  from PyQt5 import QtWidgets

  import xyzpy.utilsXyz as UXYZ
  import xyzpy.dataFromFileXyz as DFF
  import xyzpy.classFactoryXyz as CLFX
  from salomepy.onceQApplication import OnceQApplication

  verbose = False

  """
  assertEqual(a, b) 
  assertNotEqual(a, b)
  assertTrue(x)
  assertFalse(x)
  assertIs(a, b)
  assertIsNot(a, b)
  assertIsNone(x)
  assertIsNotNone(x)
  assertIn(a, b)
  assertNotIn(a, b)
  assertIsInstance(a, b)
  assertNotIsInstance(a, b)
  """

  class TestCase(unittest.TestCase):

    def test_005(self):
      ini = {'a': [1,2,3,4], 'bb': [11,22,33,44.4]}
      a = pd.DataFrame(ini)
      b = a.copy() #duplicate data
      """
      a['bb'][0] = 55.55
      generate SettingWithCopyWarning: 
      A value is trying to be set on a copy of a slice from a DataFrame
      See the caveats in the documentation: 
      http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      """
      a.loc[[0,3],'bb'] = 55.55
      a.loc[1,'bb'] = 66.66
      if verbose:
        print("\n!!!test_005\na=\n%s\nb=\n%s" % (a,b))
      self.assertTrue("55.55" in a.__str__())
      self.assertTrue("66.66" in a.__str__())
      self.assertFalse("55.55" in b.__str__())
      self.assertFalse("66.66" in b.__str__())
      self.assertTrue("11.0" in b.__str__())
      self.assertTrue("22.0" in b.__str__())
      self.assertTrue("33.0" in b.__str__())
      self.assertTrue("44.4" in b.__str__())

    def test_010(self):
      a=DFF.DataArrayXyz()
      if verbose:
        print("\n!!!test_010\n%s\n%s" % (a.__str__(), a.__repr__()))
      self.assertTrue("DataArrayXyz([])" in a.__str__())
      self.assertTrue("Empty DataFrame" in a.__repr__())
      a.setData(DFF._example_DataArray)
      if verbose:
        print("\n!!!test_010\n%s\n%s" % (a.__str__(), a.__repr__()))
      self.assertTrue("preTestScore" in a.__str__())
      self.assertTrue("preTestScore" in a.__repr__())
      self.assertFalse("Jacobson" in a.__str__())
      self.assertTrue("Jacobson" in a.__repr__())
      if verbose:
        print("\n!!!!DataArrayXyz\n")
        print(UXYZ.prettyPrintET(UXYZ.toXml(a)))
        print(UXYZ.prettyPrintET(a.toXml()))

      #test create instance from Xml
      xml = UXYZ.toXml(a)
      a2 = UXYZ.fromXml(xml)
      if verbose:
        print(UXYZ.prettyPrintET(a2.toXml()))
      self.assertEqual(a.__repr__(), a2.__repr__())


    def test_015(self):
      a=DFF.FileDataXyz("hello")
      if verbose:
        print("\n!!!test_015 %s" % a)
      self.assertEqual(a.__str__(), "hello")
      self.assertEqual(a.__repr__(), "FileDataXyz('hello')")

    def test_020(self):
      a=DFF.DataFromFileXyz()
      a.setData(DFF._example_DataArray)
      if verbose:
        print("\n!!!test_020\n%s\n%s" % (a.__str__(), a.__repr__()))
      self.assertTrue("preTestScore" in a.__repr__())
      if verbose:
        print("\n!!!!DataFromFileXyz\n")
        print(UXYZ.prettyPrintET(UXYZ.toXml(a)))
        print(UXYZ.prettyPrintET(a.toXml()))

      #test create instance from Xml
      xml = UXYZ.toXml(a)
      a2 = UXYZ.fromXml(xml)
      strXml = UXYZ.prettyPrintET(a.toXml())
      strXml2 = UXYZ.prettyPrintET(a2.toXml())
      if verbose:
        print(strXml)
        print(strXml2)
      self.assertEqual(strXml, strXml2)
      self.assertEqual(a.__repr__(), a2.__repr__())

    def test_110(self):
      import xyzpy.controllerXyz as CXYZ
      import xyzpy.guiXyz.treeXmlXyz as TXYZ
      from PyQt5 import QtGui
      app = OnceQApplication()
      C = CXYZ.ControllerXyz(desktop=QtWidgets.QMainWindow())
      V1 = TXYZ.TreeXmlXyz()
      C.setView(V1)
      M = DFF.DataFromFileXyz()
      C.setModel(M)
      if verbose:
        V1.show()
        app.exec_()
      self.assertTrue("False</DataModified>" in UXYZ.prettyPrintET(M.toXml()))


  if __name__ == '__main__':
    verbose = True
    unittest.main()


