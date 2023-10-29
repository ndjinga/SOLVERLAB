#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
import xyzpy.helpsFactoryXyz as HLFX
from xyzpy.baseXyz import _XyzConstrainBase

verbose = False

class TestCase(unittest.TestCase):
  class ExampleXyz(_XyzConstrainBase):
    _attributesList = [
      ("hello", "StrXyz"), 
      ("goodbye", "StrXyz"), 
    ]
    _helpDict = {
      "hello": ("to say hello", "hello\nbonjour..."),
      "goodbye": ("to say goodbye", "goodbye\nau revoir..."),
     }
  
  def test_010(self):
    res = HLFX.getCommonHelp("aBugNameClass", "anAttribute")
    self.assertEqual(res, None)
  
  def test_020(self):
    b = self.ExampleXyz()
    HLFX.addInCommonHelp(b)
    #print "helpsFactory Xml:\n", HLFX.toStrXml()
    res = HLFX.getCommonHelp("ExampleXyz", "hello")
    self.assertEqual(res.shortHelp, b._helpDict["hello"][0])
    self.assertEqual(res.longHelp, b._helpDict["hello"][1])
    self.assertEqual(HLFX.getCommonToolTip("ExampleXyz", "hello"), res.shortHelp)
    if verbose: print("helpsFactory Xml:\n", HLFX.toStrXml())
        
 
if __name__ == '__main__':
  verbose = False
  unittest.main()
  pass
