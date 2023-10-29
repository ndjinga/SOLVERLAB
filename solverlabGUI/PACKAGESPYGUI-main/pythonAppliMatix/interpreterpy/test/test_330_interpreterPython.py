#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import unittest
import pprint as PP

from interpreterpy.interpreterPython import InterpreterPython as IPY

verbose = False

class TestCase(unittest.TestCase):

  def test_005(self):
    ip = IPY()

    source = """\
#############################################
# an example of code source python for test #
#############################################

from __future__ import division

bb = 'hello'
cc = 111.
dd = 44
ee = cc+dd
ff = 'coucou'
_notInResult = 1 #have to be hidden
tt = 55.
ii = _notInResult/3 #have to be divide float
"""

    ip.runcode(source)
    self.assertTrue(ip.isOk())
    ipvars = ip.getVars()
    if verbose:
      print("\n########### test_005")
      print("InterpreterPython.getStdErr():\n'%s'" % ip.getStdErr())
      print("\nCode source:\n'''\n%s\n'''\n" % source)
      print("\nVariables evaluation in code:")
      for k in sorted(ipvars.keys()):
        print("%s = %s" % (k, PP.pformat(ipvars[k])))

    self.assertEqual(ipvars["bb"], 'hello')
    self.assertEqual(ipvars["cc"], 111.)
    self.assertEqual(ipvars["dd"], 44)
    self.assertEqual(ipvars["ee"], 155.)
    self.assertEqual(ipvars["ff"], 'coucou')
    self.assertFalse("_notInResult" in list(ipvars.keys()))
    self.assertEqual(ipvars["tt"], 55.)
    self.assertEqual(ipvars["ii"], 1./3.)


  def test_010(self):
    """
import code
ip = code.InteractiveInterpreter()
ip.runcode(source)
    """
    ip = IPY()
    source = """\
#############################################
# an example of code source python for test #
#############################################

from __future__ import division

dd = 'ok'
bb = there is a fait expres python syntax bug here
cc another one
"""

    """ generate StdErr error message:
File "<string>", line 8
    bb = there is a fait expres python syntax bug here
                         ^
SyntaxError: invalid syntax
    """

    ip.runcode_try(source)
    self.assertFalse(ip.isOk())
    self.assertTrue("SyntaxError: invalid syntax" in ip.getStdErr())
    if verbose:
      print("\n########### test_010")
      print("\nCode source:\n'''\n%s\n'''\n" % source)

  def test_015(self):
    ip = IPY()

    source = """\
#############################################
# an example of code source python for test #
#############################################

from __future__ import division

bb = 'hello'
dd = inexistingVar
cc = 'hello'
"""
    ip.runcode_try(source)
    self.assertFalse(ip.isOk())
    self.assertTrue("NameError: name 'inexistingVar' is not defined" in ip.getStdErr())
    ipvars = ip.getVars()
    if verbose:
      print("\n########### test_015")
      print("\nCode source:\n'''\n%s\n'''\n" % source)
      print("\nVariables evaluation in code:")
      for k in sorted(ipvars.keys()):
        print("%s = %s" % (k, PP.pformat(ipvars[k])))

    self.assertEqual(ipvars["bb"], 'hello')
    self.assertFalse("cc" in list(ipvars.keys()))



if __name__ == '__main__':
  verbose = False
  unittest.main()
  pass
