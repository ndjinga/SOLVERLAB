#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import unittest
import pprint as PP

from interpreterpy.variablesInterpreterPythonXyz import VariablesInterpreterPythonXyz as VIPY

verbose = False

class TestCase(unittest.TestCase):

  def test_005(self):
    vip = VIPY()
    if verbose: 
      print(vip.__repr__())
    self.assertEqual(vip.MyFirstVar, str(1./3.))
    source = """\
#############################################
# an example of code source python for test #
#############################################

from __future__ import division

bb = 'hello'
cc = 111.
dd = 44
ff = 1/3     #have to be divide float
ee = cc+dd
"""
    vip.setCurrentVariables(source)
    self.assertEqual(vip.bb, 'hello')
    self.assertEqual(vip.cc, str(111.))
    self.assertEqual(vip.dd, str(44))
    self.assertEqual(vip.ee, str(111.+44))
    self.assertEqual(vip.ff, str(1./3.))
    #is sorted...
    self.assertEqual(vip.getVariablesNames(), ['bb', 'cc', 'dd', 'ee', 'ff'])
    self.assertEqual(vip.getSource().split(), source.split()) #warning: last \n disappear

      

if __name__ == '__main__':
  verbose = False
  unittest.main()
  pass
