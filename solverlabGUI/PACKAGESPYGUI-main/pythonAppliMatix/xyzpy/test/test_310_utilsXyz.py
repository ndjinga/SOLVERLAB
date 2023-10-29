#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import xml.etree.ElementTree as ET

import xyzpy.utilsXyz as UXYZ
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()


verbose=False

class TestCase(unittest.TestCase):

  class Unknown(object):
    pass
  
  def test_005(self):
    a="  123   456  "
    self.assertEqual(UXYZ.stripAll(a), "123456")
    a="  123\n  456  "
    self.assertEqual(UXYZ.stripAll(a), "123\n456")
  
  def test_010(self):
    a=123
    self.assertEqual(int(UXYZ.toStrForXml(a)), 123)
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "123")
    a=123.56
    self.assertEqual(float(UXYZ.toStrForXml(a)), 123.56)
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "1.23560e+02")
    a=True
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "True") #True != .TRUE. as fortran
    self.assertEqual(bool(UXYZ.toStrForXml(a)), True)
    a=False
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "False")
    self.assertEqual(bool(UXYZ.toStrForXml(a)), True) #DANGER bool(str) True if str != ""
    a=self.Unknown()
    self.assertRaises(Exception , UXYZ.toStrForXml, a)
    
  def test_020(self):
    a=[123]
    self.assertEqual(int(UXYZ.toStrForXml(a)), 123)
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "123")
    a=[123.56]
    self.assertEqual(float(UXYZ.toStrForXml(a)), 123.56)
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "1.23560e+02")
    a=[True]
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "True") #True != .TRUE. as fortran
    self.assertEqual(bool(UXYZ.toStrForXml(a)), True)
    a=[False]
    self.assertEqual(UXYZ.toStrForXml(a).lstrip(), "False")
    self.assertEqual(bool(UXYZ.toStrForXml(a)), True) #DANGER bool(str) True if str != ""
    a=[self.Unknown()]
    self.assertRaises(Exception , UXYZ.toStrForXml, a)
    
  def test_030(self):
    nb = 12
    a=[123]*nb
    self.assertEqual(UXYZ.toStrForXml(a).split(), ["123"]*nb)
    a=[123.56]*nb
    self.assertEqual(UXYZ.toStrForXml(a).split(), ["1.23560e+02"]*nb)
    a=[True]*nb
    self.assertEqual(UXYZ.toStrForXml(a).split(), ["True"]*nb)
    a=[False]*nb
    self.assertEqual(UXYZ.toStrForXml(a).split(), ["False"]*nb) #DANGER bool(str) True if str != ""
    a=[self.Unknown()]*nb
    self.assertRaises(Exception , UXYZ.toStrForXml, a)
    
  def test_050(self):
    a={}
    verb = verbose
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_050_1:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<aDict" in res.split('\n')[1], True)
    a={"hello": 123}
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_050_2:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<hello" in res, True)
    self.assertEqual("123" in res, True)
    
    a["more"]= "456"
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_050_3:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<hello" in res, True)
    self.assertEqual("123" in res, True)
    self.assertEqual("<more" in res, True)
    self.assertEqual("456" in res, True)
    
    b={"b1": "bb1", "b2": "bb2"}
    a["b"] = b
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_050_4:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<b1" in res, True)
    self.assertEqual("bb1" in res, True)
    self.assertEqual("<b2" in res, True)
    self.assertEqual("bb2" in res, True)
    
    #bug loop
    b["loop"] = b
    self.assertRaises(Exception , UXYZ.toXml, a)
    
  def test_060(self):
    a=[]
    verb = verbose
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_060_1:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<aList" in res.split('\n')[1], True)
    a.append(123)
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_060_2:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<_0_" in res, True)
    self.assertEqual("123" in res, True)
    
    a.append("456")
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_060_3:\n%s" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<_0_" in res, True)
    self.assertEqual("123" in res, True)
    self.assertEqual("<_1_" in res, True)
    self.assertEqual("456" in res, True)
    
    b={"b1": "bb1", "b2": "bb2"}
    a.append(b)
    res = UXYZ.prettyPrintET(UXYZ.toXml(a))
    if verb: logger.info("test_060_4:\n" % UXYZ.prettyPrintET(UXYZ.toXml(a)))
    self.assertEqual("<b1" in res, True)
    self.assertEqual("bb1" in res, True)
    self.assertEqual("<b2" in res, True)
    self.assertEqual("bb2" in res, True)
    
    """
test_060_4:
<?xml version='1.0' encoding='UTF-8'?>
<aList size='3' typeClass='list'>
  <_0_ typeClass='IntXyz' index='0'>          123</_0_>
  <_1_ typeClass='StrXyz' index='1'>          456</_1_>
  <_2_ typeClass='dict' index='2'>
    <b1 typeClass='StrXyz'>          bb1</b1>
    <b2 typeClass='StrXyz'>          bb2</b2>
  </_2_>
</aList>
"""
     
    #bug loop
    b["loop"] = b
    self.assertRaises(Exception , UXYZ.toXml, a)
    
  def test_070(self):
    strXml = """<?xml version='1.0' encoding='UTF-8'?>
<aList size='3' typeClass='list'>
  <_0_ typeClass='IntXyz' index='0'>          123</_0_>
  <_1_ typeClass='StrXyz' index='1'>          456</_1_>
  <_2_ typeClass='dict' index='2'>
    <b1 typeClass='StrXyz'>          bb1</b1>
    <b2 typeClass='StrXyz'>          bb2</b2>
  </_2_>
</aList>"""
    aTruc = UXYZ.fromXml(strXml)
    if verbose: logger.info("aTruc\n%s" % aTruc.__repr__())
    self.assertEqual(aTruc[0], 123)
    self.assertEqual(aTruc[1], "456")
    self.assertEqual(aTruc[2]["b1"], "bb1")
    self.assertEqual(aTruc[2]["b2"], "bb2")

  def test_100(self):
    root = ET.Element('Root')
    res = UXYZ.prettyPrintET(root)
    if verbose: logger.info("test_040_0\n%s" % res)
    self.assertEqual(res[0:14], "<?xml version=")
    self.assertEqual(res.split('\n')[1], "<Root/>")
    
    root.text = UXYZ.toStrForXml(123)
    res = UXYZ.prettyPrintET(root)
    if verbose: logger.info("test_040_1\n%s" % res)
    self.assertEqual(UXYZ.stripAll(res.split('\n')[1]), "<Root>123</Root>")
    
    root.text = UXYZ.toStrForXml([123]*25)
    res = UXYZ.prettyPrintET(root)  
    if verbose: logger.info("test_040_2\n%s" % res)
    self.assertEqual(len(res.split(' 123')), 26)
   
   
if __name__ == '__main__':
  unittest.main()
  pass
