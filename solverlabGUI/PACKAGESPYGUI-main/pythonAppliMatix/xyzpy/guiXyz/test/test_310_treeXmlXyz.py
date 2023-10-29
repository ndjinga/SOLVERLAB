#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import sys
import unittest
import pprint as PP

from PyQt5 import QtCore, QtGui, QtWidgets
import xml.etree.ElementTree as ET
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem, TreeXmlXyzMainWidget
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX


setTimer = True
deltaTime = 500
withShow = True
verbose = False
user = os.getenv("USERNAME")
userx = ""


class TestCase(unittest.TestCase):
  xmlStr = """\
<?xml version="1.0" encoding="utf-8"?>
<data>
  <country name="Liechtenstein">
    <rank>1</rank>
    <year>2008</year>
    <etc>141100</etc>
    <neighbor direction="E" name="Austria"/>
    <neighbor direction="W" name="Switzerland"/>
  </country>
  <country name="Singapore">
    <rank>3</rank>
    <year>2011</year>
    <etc>59900</etc>
    <neighbor direction="N" name="Malaysia"/>
  </country>
  <country name="Panama">
    <rank>68</rank>
    <year>2011</year>
    <etc>13600</etc>
    <neighbor direction="W" name="Costa Rica"/>
    <neighbor direction="E" name="Colombia"/>
  </country>
</data>
"""
  old_xmlStr = """\
<?xml version="1.0" encoding="utf-8"?>
<data>
  <country name="Liechtenstein">
    <rank typeClass="IntXyz">1</rank>
    <year typeClass="FloatPosXyz">2008</year>
    <gdppc typeClass="FloatXyz">141100</gdppc>
    <neighbor direction="E" name="Austria"/>
    <neighbor direction="W" name="Switzerland"/>
  </country>
  <country name="Singapore">
    <rank typeClass="Int04Xyz">3</rank>
    <year typeClass="IntPosXyz">2011</year>
    <gdppc typeClass="FloatXyz">59900</gdppc>
    <neighbor direction="N" name="Malaysia"/>
  </country>
  <country name="Panama">
    <rank>68</rank>
    <year>2011</year>
    <gdppc typeClass="FloatSupEq1Xyz">13600</gdppc>
    <neighbor direction="W" name="Costa Rica"/>
    <neighbor direction="E" name="Colombia"/>
  </country>
</data>
"""
  def root_path():
    return os.path.abspath(os.sep)
  
  def new_element(self, name, text, kwargs=None):
    """kwargs are as xml tag attribute"""
    tag = ET.Element(name)
    #tag.text = UXYZ.toStrForXml(text) #just a little more human readable
    tag.text = text #just almost human readable
    if kwargs != None: tag.attrib = kwargs
    return tag
  
  def launchTimer(self, app, wid):
    timer = QtCore.QTimer();
    timer.timeout.connect(wid.close)
    if setTimer: timer.start(deltaTime)
    app.exec_()
    
  def test_010(self):
    a = TreeXmlXyzItem()
    self.assertEqual(a.filterToStr(None), ("", ""))
    self.assertEqual(a.filterToStr("essai"), ("essai", "essai"))
    self.assertEqual(a.filterToStr("essai\n"), ("essai", "essai"))
    c40 = "0000000000111111111122222222223333333333"
    c30 = "000000000011111111112222222222"
    c30p = "000000000011111111112222222222..."
    self.assertEqual(a.filterToStr(c40), (c30p, c40))
    self.assertEqual(a.filterToStr(c30), (c30, c30))
    self.assertEqual(a.filterToStr("\n"), ("", ""))
    self.assertEqual(a.filterToStr("\n"*40), ("", ""))
    lines2 = " line1\n line2\n"
    lines3 = " line1\n line2\n line3\n"
    self.assertEqual(a.filterToStr(lines2), ('line1...', 'line1\n line2'))
    self.assertEqual(a.filterToStr(lines3), ('line1...', 'line1\n line2\n line3'))
    d2 = {"k1": "ak1key", "k2": "ak2key"}
    d3 = {"k1": "ak1key", "k2": "ak2key", "k3": "ak3key"}
    self.assertEqual(a.filterToStr(d2), ("k1='ak1key'...", "k1='ak1key'\nk2='ak2key'"))
    self.assertEqual(a.filterToStr(d3), ("k1='ak1key'...", "k1='ak1key'\nk2='ak2key'\nk3='ak3key'"))
    if verbose: #verbose test ,not for allTestLauncher.py
      self.assertEqual("?type text?" in a.filterToStr(a)[0], True)

  def test_100(self):
    if userx in ["wambeke", "matix"]:
      print("\nWARNING: test_100 skipped as %s" % user)
      return
    app = OnceQApplication()
    data = ET.fromstring(self.xmlStr)
    treeWid = TreeXmlXyz()
    treeWid.setFromXml(data)
    treeWid.resize(500, 500)
    treeWid.show()
    self.launchTimer(app, treeWid)
  
  def test_110(self):
    if userx in ["wambeke", "matix"]:
      print("\nWARNING: test_110 skipped as %s" %user)
      return
    app = OnceQApplication()
    data = ET.fromstring(self.xmlStr)
    treeWid = TreeXmlXyzMainWidget()
    treeWid.widget.setFromXml(data)
    treeWid.resize(500, 500)
    treeWid.show()
    self.launchTimer(app, treeWid)
    
  def test_200(self):
    if userx in ["wambeke", "matix"]:
      print("\nWARNING: test_200 skipped as %s" %user)
      return
    app = OnceQApplication()
    allClasses = CLFX.getAllXyzClasses()
    tag = ET.Element('allClasses')
    res = []
    ok = True
    for aNameClass in sorted(allClasses.keys()):
      if verbose: print("test_200 nameClass %s" % aNameClass)
      aClass = allClasses[aNameClass]
      anInstance = aClass()  # default init
      try:
        # only for simple immutables
        AStrValue = UXYZ.toStrForXml(anInstance)
        tag.append(self.new_element(aNameClass, AStrValue, {"typeClass": aNameClass}))
      except Exception as e:
        # no direct representation of a BaseFreeXyz for example
        # AStrValue = anInstance.__repr__()
        AStrValue = str(e)
        if "works only for immutables" not in AStrValue:
          tag.append(self.new_element(aNameClass, AStrValue, {"typeClass": aNameClass}))
          ok = False

      res.append([aNameClass, AStrValue])

    if not ok: print("test_200\n%s" % PP.pformat(res, width=150))
    treeWid = TreeXmlXyzMainWidget()
    treeWid.widget.setFromXml(tag)
    treeWid.resize(600, 400)
    treeWid.show()
    if verbose:
      app.exec()
    else:
      self.launchTimer(app, treeWid)
    self.assertTrue(ok)

if __name__ == '__main__':
  verbose = True
  setTimer = False #user call and wait
  user=""
  #from xyzpy.guiXyz import treeXmlXyz
  #treeXmlXyz.verboseEvent = not setTimer
  unittest.main()
  pass
