#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import sys
import glob
import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
import xml.etree.ElementTree as ET
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.dialogXmlXyz import DialogXmlXyz
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX


setTimer = True
deltaTime = 2000
withShow = True
verbose = False
user = os.getenv("USERNAME")

# reading file from children directory ./tests_ekinox
PACKAGESPY_ROOT_DIR = os.getenv("PACKAGESPY_ROOT_DIR")
if PACKAGESPY_ROOT_DIR is not None:
  testDir = os.path.realpath( os.path.join(PACKAGESPY_ROOT_DIR, "pythonAppliMatix", "ekinoxpy", "tests_ekinox") )
  if verbose: print("testDir %s" % testDir)


class TestCase(unittest.TestCase):
  #it is only for test example, hand modified, it is NOT EkinoxEki example
  xmlStr = """\
<?xml version="1.0" encoding="utf-8"?>
<EkinoxEki treePyName="" typeClass="BaseFreeXyz">
  <release treePyName=".release" typeClass="StrXyz">          1.0</release>
  <AlliageEki treePyName=".AlliageEki" typeClass="BaseFreeXyz">
    <general treePyName=".AlliageEki.general" typeClass="BaseFreeXyz">
      <launchMode treePyName=".AlliageEki.general.launchMode" typeClass="StrXyz">      initial</launchMode>
      <restartFile treePyName=".AlliageEki.general.restartFile" typeClass="FileXyz">             </restartFile>
      <boxSize treePyName=".AlliageEki.general.boxSize" typeClass="FloatPosXyz">  2.00000e+02</boxSize>
      <numberOfLayersInMetal treePyName=".AlliageEki.general.numberOfLayersInMetal" typeClass="IntPosXyz">          200</numberOfLayersInMetal>
      <numberOfLayersInOxyd treePyName=".AlliageEki.general.numberOfLayersInOxyd" typeClass="IntPosXyz">            3</numberOfLayersInOxyd>
    </general>
    <metal treePyName=".AlliageEki.metal" typeClass="BaseFreeXyz">
      <name treePyName=".AlliageEki.metal.name" typeClass="StrXyz">       FeNiCr</name>
      <stoechiometry treePyName=".AlliageEki.metal.stoechiometry" typeClass="StrXyz">        Cr2O3</stoechiometry>
      <CNi treePyName=".AlliageEki.metal.CNi" typeClass="Float01Xyz">  0.00000e+00</CNi>
      <CCr treePyName=".AlliageEki.metal.CCr" typeClass="Float01Xyz">  2.20000e-01</CCr>
      <CFe treePyName=".AlliageEki.metal.CFe" typeClass="Float01Xyz">  7.80000e-01</CFe>
      <rho treePyName=".AlliageEki.metal.rho" typeClass="BaseFreeXyz">
        <rhoVolume treePyName=".AlliageEki.metal.rho.rhoVolume" typeClass="FloatPosXyz">  1.00000e+04</rhoVolume>
        <rhoSurface treePyName=".AlliageEki.metal.rho.rhoSurface" typeClass="FloatPosXyz">  1.00000e+10</rhoSurface>
        <rhoInterface treePyName=".AlliageEki.metal.rho.rhoInterface" typeClass="FloatPosXyz">  1.00000e+12</rhoInterface>
      </rho>
      <gamma treePyName=".AlliageEki.metal.gamma" typeClass="FloatPosXyz">  1.50000e+00</gamma>
      <pbr treePyName=".AlliageEki.metal.pbr" typeClass="FloatPosXyz">  2.00000e+00</pbr>
    </metal>
    <oxyd treePyName=".AlliageEki.oxyd" typeClass="BaseFreeXyz">
      <profil treePyName=".AlliageEki.oxyd.profil" typeClass="StrXyz">         Flat</profil>
      <anionicConcentrations treePyName=".AlliageEki.oxyd.anionicConcentrations" typeClass="BaseFreeXyz">
        <metalOxyd treePyName=".AlliageEki.oxyd.anionicConcentrations.metalOxyd" typeClass="FloatPosXyz">  1.00000e-04</metalOxyd>
        <oxydGaz treePyName=".AlliageEki.oxyd.anionicConcentrations.oxydGaz" typeClass="FloatPosXyz">  7.16000e-04</oxydGaz>
        <intermediate treePyName=".AlliageEki.oxyd.anionicConcentrations.intermediate" typeClass="FloatPosXyz">  5.00000e-04</intermediate>
      </anionicConcentrations>
      <cationicConcentrations treePyName=".AlliageEki.oxyd.cationicConcentrations" typeClass="BaseFreeXyz">
        <metalOxyd treePyName=".AlliageEki.oxyd.cationicConcentrations.metalOxyd" typeClass="FloatPosXyz">  1.00000e-06</metalOxyd>
        <oxydGaz treePyName=".AlliageEki.oxyd.cationicConcentrations.oxydGaz" typeClass="FloatPosXyz">  1.00000e-06</oxydGaz>
        <intermediate treePyName=".AlliageEki.oxyd.cationicConcentrations.intermediate" typeClass="FloatPosXyz">  5.00000e-03</intermediate>
      </cationicConcentrations>
      <alphanum treePyName=".AlliageEki.oxyd.alphanum" typeClass="FloatPosXyz">  0.00000e+00</alphanum>
    </oxyd>
    <oxydationParameters treePyName=".AlliageEki.oxydationParameters" typeClass="BaseFreeXyz">
      <temperature treePyName=".AlliageEki.oxydationParameters.temperature" typeClass="FloatPosXyz">  1.53300e+03</temperature>
      <timeStep treePyName=".AlliageEki.oxydationParameters.timeStep" typeClass="FloatPosXyz">  4.00000e-01</timeStep>
      <nti treePyName=".AlliageEki.oxydationParameters.nti" typeClass="IntPosXyz">     10000000</nti>
    </oxydationParameters>
  </AlliageEki>
</EkinoxEki>
"""

  def root_path():
    return os.path.abspath(os.sep)

  def new_element(self, name, text, kwargs=None):
    """kwargs are as xml tag attribute"""
    tag =  ET.Element(name)
    #tag.text = UXYZ.toStrForXml(text) #just a little more human readable
    tag.text = text #just almost human readable
    if kwargs != None: tag.attrib = kwargs
    return tag

  def launchTimer(self, app, wid):
    timer = QtCore.QTimer();
    timer.timeout.connect(wid.close)
    if setTimer: timer.start(deltaTime)
    app.exec_()

  def test_100(self):
    app = OnceQApplication()
    data = None
    wid = DialogXmlXyz()
    wid.setFromXml(data)
    #wid.activeStretch()
    wid.resize(500, 500)
    wid.show()
    self.launchTimer(app, wid)

  def test_200(self):
    app = OnceQApplication()
    data = ET.fromstring(self.xmlStr)
    localModel = UXYZ.fromXml(data)
    alliageData = localModel.AlliageEki
    data = alliageData.toStrXml()
    if verbose: print("test dialog on: \n",data)
    wid = DialogXmlXyz()
    wid.setFromXml(data)
    #wid.activeStretch()
    wid.resize(500, 500)
    wid.show()
    self.launchTimer(app, wid)


if __name__ == '__main__':
  setTimer = False #user call and wait
  unittest.main()
  pass
