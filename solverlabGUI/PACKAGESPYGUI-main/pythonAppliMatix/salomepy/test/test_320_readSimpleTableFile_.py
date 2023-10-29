#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import unittest
import salomepy.readSimpleTableFile as RSTF
import salomepy.utilsWorkdir as UTW
runningDir = UTW.getWorkdirDefault("TESTS")  

verbose = False
_aName = "simpleTableFileTest.csv"
_aFile = os.path.join(runningDir, _aName)

class TestCase(unittest.TestCase):
  
  def test_010(self):
    s = ""
    for i in range(10):
      ii = float(i)
      s += "%i %.5e %.4f %s\n" % (ii*2+1, (ii+10)/2, ii**1.5, str((11-ii)**.5))
    with open(_aFile, "w") as F: F.write(s)
    xy, titles_xy, aName = RSTF.readSimpleTableFile( _aFile )
    self.assertEqual(len(titles_xy), 4)
    self.assertEqual(len(xy), 4)
    self.assertEqual(len(xy[0]), 10)
    self.assertEqual(aName, _aName)

  def test_020(self):
    s = """\
#aa bb cc
11 22 33 44 #only 3 columns on header
55 66 77 88
"""
    with open(_aFile, "w") as F: F.write(s)
    xy, titles_xy, aName = RSTF.readSimpleTableFile( _aFile )
    self.assertEqual(len(titles_xy), 3)
    self.assertEqual(titles_xy, ['aa','bb','cc'])
    self.assertEqual(len(xy), 3)
    self.assertEqual(len(xy[0]), 2)
    self.assertEqual(xy[0], [11., 55.])
    self.assertEqual(aName, _aName)

  def test_030(self):
    s = """\
#  aa    bb cc dd  
11 22 33e0 44
55 66 77 88D2 #double fortran
"""
    with open(_aFile, "w") as F: F.write(s)
    xy, titles_xy, aName = RSTF.readSimpleTableFile( _aFile )
    self.assertEqual(len(titles_xy), 4)
    self.assertEqual(titles_xy, ['aa','bb','cc','dd'])
    self.assertEqual(len(xy), 4)
    self.assertEqual(len(xy[0]), 2)
    self.assertEqual(xy[0], [11., 55.])
    self.assertEqual(xy[3], [44., 8800.])
    self.assertEqual(aName, _aName)

    
if __name__ == '__main__':
  unittest.main()
  


