#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
from salomepy.onceQApplication import OnceQApplication
from PyQt5.QtGui import QIcon
import salomepy.iconsUser as IUSR
import xyzpy.loggingXyz as LOG

class TestCase(unittest.TestCase):
  def test_000(self):
    #need an QApplication
    app=OnceQApplication() #warning:  problems for ParexIcons with ubuntu
    noIcon = IUSR.getIconFromName("noIcon")
    self.assertNotEqual(noIcon, None)
    self.assertEqual(noIcon.__class__, QIcon().__class__)
    LOG.pushLevel("CRITICAL")
    self.assertEqual(IUSR.getIconFromName("inexistingToto"), noIcon)
    LOG.popLevel()

  def test_010(self):
    #need an QApplication
    app=OnceQApplication() #warning:  problems for ParexIcons with ubuntu
    noIcon = IUSR.getIconFromName("noIcon")
    for i in list(IUSR.IconsUser.keys()):
      icon = IUSR.getIconFromName(i)
      if i == "noIcon":
        self.assertEqual(icon, noIcon)

  def test_020(self):
    #need an QApplication
    app=OnceQApplication() #warning:  problems for ParexIcons with ubuntu
    noIcon = IUSR.getIconFromName("noIcon")
    icon = IUSR.getIconFromPattern("inexistingPattern")
    self.assertNotEqual(icon, noIcon)
    self.assertEqual(icon, None)
    for pat, name in IUSR.IconsUser_fnpatterns:
      icon = IUSR.getIconFromPattern(pat)
      self.assertNotEqual(icon, noIcon)
      self.assertNotEqual(icon, None)

if __name__ == '__main__':
  """IconsUserTest.py"""
  unittest.main()
  pass
