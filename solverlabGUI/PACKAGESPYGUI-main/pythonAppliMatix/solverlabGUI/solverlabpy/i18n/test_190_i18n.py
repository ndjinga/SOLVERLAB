#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2018  CEA/DEN
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END

import os
import gettext
import unittest

verbose = False

class TestCase(unittest.TestCase):
  
  def test_001(self):
    # first load resources for internationalization
    gettext.install('solverlabGUI', os.path.realpath(os.path.dirname(__file__)))
 
  def test_005(self):
    res = _("Harvey writes '%(1)s' for %(2)s.") % {"1": "hello", "2": "test_005"}
    if verbose: print(res)
    self.assertEqual(res, "pour test_005 Hervé écrit 'hello'.")

  def test_010(self):
    res = _("Harvey writes '%(1)s' for %(2)s.") % {"1": _("hello"), "2": "test_010"}
    if verbose: print(res)
    self.assertEqual(res, "pour test_010 Hervé écrit 'bonjour'.")

  def test_020(self):
    # keep Ooops inexisting in solverlabGui.po as no translation
    res = _("Harvey writes '%(1)s' for %(2)s.") % {"1": _("Ooops"), "2": "test_020"}
    if verbose: print(res)
    self.assertEqual(res, "pour test_020 Hervé écrit 'Ooops'.")

if __name__ == '__main__':
  verbose = False
  unittest.main()
  pass
