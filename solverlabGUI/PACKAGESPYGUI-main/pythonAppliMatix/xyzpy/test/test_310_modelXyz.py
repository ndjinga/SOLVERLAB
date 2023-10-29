#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import glob
import unittest
import xyzpy.test.modelSimpleTest as MXYZ
import xml.etree.ElementTree as ET
import xyzpy.utilsXyz as UXYZ
import xyzpy.intFloatListXyz as ILFX

verbose=False

class TestCase(unittest.TestCase):
  def test_010(self):
    M = MXYZ.ModelSimpleTest()
    self.assertEqual("ModelSimpleTest[" in M.objectName(), True)
    self.assertEqual(M.lastReceiveRequest(), None)
    self.assertRaises(Exception , M.receiveRequest, "n'importe quoi!" )
    strXmlRequest = """<?xml version='1.0' encoding='UTF-8'?>
    <BaseFreeXyz typeClass='BaseFreeXyz'>
      <typeRequest typeClass='StrXyz'>hello...</typeRequest>
    </BaseFreeXyz>
    """
    M.receiveRequest(strXmlRequest)
    self.assertEqual(M.lastReceiveRequest(), strXmlRequest)
    strXmlRequestWithBug = """<?xml version='1.0' encoding='UTF-8'?>
    <BaseFreeXyz typeClass='BaseFreeXyz'>
      <typeRequest typeClass='StrXyz'>hello...</typeNASTYBUGRequest>
    </BaseFreeXyz>
    """
    self.assertRaises(Exception , M.receiveRequest, strXmlRequestWithBug)

if __name__ == '__main__':
  unittest.main()
  pass
