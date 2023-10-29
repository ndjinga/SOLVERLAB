#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import glob
import unittest
import xyzpy.baseXyz as BXYZ
import xml.etree.ElementTree as ET
import xyzpy.utilsXyz as UXYZ
import xyzpy.intFloatListXyz as IFLX

verbose=False

########################################################################################
class TestCase(unittest.TestCase):
  
  def test_005(self):
    a = BXYZ.BaseXyz()
    a.tintin = "reporter"
    a.milou = "dog"
    self.assertEqual(a.tintin, "reporter")
    self.assertEqual(a.milou, "dog")
    del(a.tintin)
    self.assertFalse(hasattr(a, "tintin"))
    self.assertEqual(a.milou, "dog")

  def test_006(self):
    a = []
    self.assertEqual(len(a), 0)
    a.append(1)
    self.assertEqual(len(a), 1)
    a.append(2)
    self.assertEqual(len(a), 2)
    self.assertEqual(a, [1,2])
    del(a[0])
    self.assertEqual(a, [2])
    del(a[0])
    self.assertEqual(a, [])
    a=[1,2,3]
    a.insert(1,44)
    self.assertEqual(a, [1, 44, 2, 3])
    a.insert(10,55)
    self.assertEqual(a, [1, 44, 2, 3, 55])
    del(a[1])
    self.assertEqual(a, [1, 2, 3, 55])
    self.assertEqual(len(a), 4)
    
  def test_007(self):
    a = BXYZ.ListOfBaseXyz()
    self.assertEqual(len(a), 0)
    self.assertRaises(Exception , a.append, 1 )
    a.append(IFLX.IntXyz(1))
    self.assertEqual(len(a), 1)
    a.append(IFLX.IntXyz(2))
    self.assertEqual(len(a), 2)
    self.assertEqual(a, [1,2])
    del(a[0])
    self.assertEqual(a, [2])
    del(a[0])
    self.assertEqual(a, [])
    a=BXYZ.ListOfBaseXyz()
    a.append(IFLX.IntXyz(1))
    a.append(IFLX.IntXyz(2))
    a.append(IFLX.IntXyz(3))
    a.insert(1,IFLX.IntXyz(44))
    self.assertEqual(a, [1, 44, 2, 3])
    a.insert(10,IFLX.IntXyz(55))
    self.assertEqual(a, [1, 44, 2, 3, 55])
    del(a[1])
    self.assertEqual(a, [1, 2, 3, 55])
    self.assertEqual(len(a), 4)

  def test_010(self):
    a = BXYZ.BaseFreeXyz()
    self.assertEqual(a.isSet(), True)
    self.assertRaises(Exception , a.setIsSet, "bug" )
    a.setIsSet(False) #useless for Free, useful for inherited
    self.assertEqual(a.isSet(), False)
    a.setIsSet(True)
    self.assertEqual(a.isSet(), True)
    
    self.assertEqual(a.getAttributes(), [])
    self.assertEqual(a.getCurrentAttributes(), [])
    self.assertRaises(Exception , a.__setattr__, "toto", 1 )
    
  def test_020(self):
    a = BXYZ.BaseFreeXyz()
    b = BXYZ.BaseFreeXyz()
    
    #it is a choice important do not change it
    #a == b only if ids are equals
    self.assertNotEqual(a, b) 
    self.assertEqual(a.equal(b), True)
    self.assertEqual(a.equal(None), False)
    self.assertRaises(Exception , a.equal,  1)
    
  def test_025(self):
    a = BXYZ.BaseFreeXyz()
    a1 = IFLX.IntXyz(1)
    a2 = IFLX.IntXyz(2)
    a3 = IFLX.IntXyz(3)
    a.att1 = a1
    a.att2 = a2
    ida2 = id(a2)
    
    self.assertEqual(a.att1, 1)
    self.assertEqual(a.att2, 2)
    res = [ai for ai in a]
    self.assertEqual(res, [a1,  a2]) 
    res = [id(ai) for ai in a]
    self.assertEqual(res, [id(a1),  id(a2)]) 
    
    self.assertEqual(a2.parentAsAttribute(), a) 
    #delattr(a, 'att2')
    a.__delattr__('att2')
    self.assertEqual(id(a2), ida2) 
    self.assertEqual(a2.parentAsAttribute(), None) 
    
    res = [ai for ai in a]
    self.assertEqual(res, [a1]) 
    a.att2 = a3
    self.assertEqual(a3.parentAsAttribute(), a) 
    self.assertEqual(a2, 2) 
    self.assertEqual(a.att2, 3) 
    res = [ai for ai in a]
    self.assertEqual(res, [a1,  a3]) 
    
    a.att2 = a2
    self.assertEqual(a3.parentAsAttribute(), None) 
    self.assertEqual(a3, 3) 
    self.assertEqual(a2, 2) 
    self.assertEqual(a.att2, 2) 
    res = [ai for ai in a]
    self.assertEqual(res, [a1,  a2]) 
    
    self.assertRaises(Exception , a.equal,  1)
    
  def test_030(self):
    a = BXYZ.BaseFreeXyz()
    b = BXYZ.BaseFreeXyz()
    
    a.toto = b
    self.assertEqual("toto" in list(a.__dict__.keys()), True)
    self.assertEqual(a.getCurrentAttributes(), ['toto'])
    self.assertEqual(a.getAttributes(), ['toto'])
        
    self.assertEqual(id(a.toto), id(b))
    self.assertEqual(id(a), id(b.parentAsAttribute()))
    self.assertRaises(Exception , a.__setattr__, "titi", b ) #avoid loop in tree

  def test_040(self):
    a = BXYZ.BaseFreeXyz()
    b = BXYZ.BaseFreeXyz()
    c = BXYZ.BaseFreeXyz()
    a.toto = b
    a.toto.tutu = c
    
    self.assertEqual("tutu" in list(b.__dict__.keys()), True)
    self.assertEqual(b.getCurrentAttributes(), ['tutu'])
    self.assertEqual(b.getAttributes(), ['tutu'])
        
    self.assertEqual(id(a.toto.tutu), id(c))
    self.assertEqual(id(a), id(c.parentAsAttribute().parentAsAttribute()))
    self.assertEqual(None, a.parentAsAttribute())
    self.assertEqual(a, a.getRoot())
    self.assertEqual(a, b.getRoot())
    self.assertEqual(a, c.getRoot())
    self.assertEqual(a, a.getRoot())
    self.assertEqual(a, a.toto.getRoot())
    self.assertEqual(a, a.toto.tutu.getRoot())
    
    res = [name for name, _  in list(a.items())]
    self.assertEqual(res , ["toto"])
    
    res = [val for val in a]
    self.assertEqual(res , [b])
    
    res = [name for name, _  in list(a.toto.items())]
    self.assertEqual(res , ["tutu"])
    
    res = [val for val in a.toto]
    self.assertEqual(res , [c])

  def test_050(self):
    a = BXYZ.BaseFreeXyz()
    b = BXYZ.BaseFreeXyz()
    c = BXYZ.BaseFreeXyz()
    a.toto = b
    a.toto.tutu = c
    
    self.assertEqual(a.toto.tutu.getNameAsAttribute() , "tutu")
    self.assertEqual(a.toto.getNameAsAttribute() , "toto")
    self.assertEqual(a.toto.tutu.getTreePyName() , ".toto.tutu")
    self.assertEqual(a.getTreePyName() , "")
    
    root = BXYZ.BaseFreeXyz()
    root.titi = a
    self.assertEqual(a.toto.tutu.getNameAsAttribute() , "tutu")
    self.assertEqual(a.toto.getNameAsAttribute() , "toto")
    self.assertEqual(a.toto.tutu.getTreePyName() , ".titi.toto.tutu")
    self.assertEqual(a.getTreePyName() , ".titi")
    self.assertEqual(root.getTreePyName() , "")
    
    #light python meta programmation
    res = {"root": root}
    aMetaStr = 'aTutu1=root.titi.toto.tutu'
    exec(aMetaStr, res)
    self.assertEqual(res["aTutu1"] , c)

    aMetaStr = 'aTutu2=root' + c.getTreePyName()
    exec(aMetaStr, res)
    self.assertEqual(res["aTutu2"] , c)
    self.assertNotEqual(res["aTutu2"] , root)

  def test_070(self):
    a = BXYZ.BaseFreeXyz()
    a.aa = IFLX.IntXyz(1)
    a.bb = IFLX.IntXyz(2)
    c = BXYZ.BaseFreeXyz()
    c.aa = IFLX.IntXyz(3)
    c.bb = IFLX.IntXyz(4)
    a.cc = c
    #print "getDictOfTreePyNameStrValues",a.getDictOfTreePyNameStrValues()
    aDict = a.getDictOfTreePyNameStrValues()
    self.assertEqual(aDict['.aa'], '1')
    self.assertEqual(aDict['.bb'], '2')
    self.assertEqual(aDict['.cc.aa'], '3')
    self.assertEqual(aDict['.cc.bb'], '4')
    
  def test_075(self):
    a = BXYZ.ListOfBaseXyz()
    a.append(IFLX.IntXyz(1))
    a.append(IFLX.IntXyz(2))
    self.assertEqual(a[0], 1)
    self.assertEqual(a[1], 2)
    a.pop()
    self.assertEqual(len(a), 1)
    a.append(IFLX.IntXyz(3))
    self.assertEqual(len(a), 2)
    self.assertEqual(a[1], 3)
    a[1] = a[1].sameType("1234")
    self.assertEqual(a[1], 1234)
    #print a,a.getDictOfTreePyNameStrValues()
    
  def test_078(self):
    
    #for test
    class ListOfAllowedBaseXyz_1(BXYZ.ListOfBaseXyz):
      _allowedClasses = [IFLX.IntXyz]
      pass
      
    class ListOfAllowedBaseXyz_2(BXYZ.ListOfBaseXyz):
      _allowedClasses = [IFLX.IntXyz, ListOfAllowedBaseXyz_1]
      pass
    
    class ListOfAllowedBaseXyz_3(BXYZ.ListOfBaseXyz):
      _allowedClasses = [IFLX.IntXyz, BXYZ.ListOfBaseXyz]
      pass
    
    a = ListOfAllowedBaseXyz_1()
    a.append(IFLX.IntXyz(1))
    self.assertEqual(a[0], 1)
    
    #a.append(IFLX.FloatXyz(1))
    x = IFLX.FloatXyz(1)
    self.assertRaises(Exception, a.append, x )
    
    b = ListOfAllowedBaseXyz_2()
    b.append(IFLX.IntXyz(1))
    self.assertEqual(b[0], 1)
    x = IFLX.FloatXyz(1)
    self.assertRaises(Exception, b.append, x)
    b.append(a)
    self.assertEqual(b[1], a)
    
    c = ListOfAllowedBaseXyz_3()
    c.append(IFLX.IntXyz(1))
    self.assertEqual(c[0], 1)
    x = IFLX.FloatXyz(1)
    self.assertRaises(Exception, c.append, x)
    x = BXYZ.ListOfBaseXyz()
    c.append(x)
    self.assertEqual(c[1], x)
    
  def test_080(self):
    a = BXYZ.ListOfBaseXyz()
    a.append(IFLX.IntXyz(1))
    a.append(IFLX.IntXyz(2))
    #print "getDictOfTreePyNameStrValues",a.getDictOfTreePyNameStrValues()
    aDict = a.getDictOfTreePyNameStrValues()
    self.assertEqual(aDict['[0]'], '1')
    self.assertEqual(aDict['[1]'], '2')

    a = BXYZ.BaseFreeXyz()
    a.aa = IFLX.IntXyz(0)
    d = BXYZ.ListOfBaseXyz()
    d.append(IFLX.IntXyz(1))
    d.append(IFLX.IntXyz(2))
    a.dd = d
    aDict = a.getDictOfTreePyNameStrValues()
    self.assertEqual(aDict['.aa'], '0')
    self.assertEqual(aDict['.dd[0]'], '1')
    self.assertEqual(aDict['.dd[1]'], '2')
    a.aa = IFLX.IntXyz(3)
    aDict = a.getDictOfTreePyNameStrValues()
    self.assertEqual(aDict['.aa'], '3')
    self.assertEqual(aDict['.dd[1]'], '2')
    
    aDict['.aa']= '33'
    aDict['.dd[1]']= '22'
    a.setValueFromDictOfTreePyNameStrValues(aDict)
    aDict = a.getDictOfTreePyNameStrValues()
    self.assertEqual(aDict['.aa'], '33')
    self.assertEqual(aDict['.dd[1]'], '22')
    
  def test_100(self):
    a = BXYZ.BaseFreeXyz()
    self.assertEqual('<BaseFreeXyz' in a.toStrXml(), True)
    b = BXYZ.BaseFreeXyz()
    c = BXYZ.BaseFreeXyz()
    a.toto = b
    a.toto.tutu = c
    aStrXml = a.toStrXml()
    #if verbose: print "test_100_0\n", aStrXml 
    self.assertEqual('<BaseFreeXyz' in aStrXml, True)
    self.assertEqual('<toto' in aStrXml, True)
    self.assertEqual('<tutu' in aStrXml, True)

  def test_200(self):
    a = BXYZ.ListOfBaseXyz() 
    self.assertEqual(len(a), 0)
    self.assertEqual(a.getNameObject(), "ListOfBaseXyz") #ListOfUnknownName")
    
    self.assertEqual(a.isSet(), True)
    self.assertRaises(Exception , a.setIsSet, "bug" )
    a.setIsSet(False) #useless here, useful for inherited
    self.assertEqual(a.isSet(), False)
    a.setIsSet(True)
    self.assertEqual(a.isSet(), True)
    
    self.assertRaises(Exception , a.__setattr__, "toto", 1 )
    
  def test_211(self):
    a = BXYZ.ListOfBaseXyz()
    b = BXYZ.BaseFreeXyz()
    a.append(b)
    self.assertEqual(len(a), 1)  
    self.assertRaises(Exception , a.append, b ) #b is appended yet somewhere
    self.assertEqual(a[0], b)
    self.assertEqual(a.getNameObject(), "ListOfBaseFreeXyz")
    c = a.pop()
    self.assertEqual(len(a), 0)
    self.assertEqual(c, b)
    
  def test_210(self):
    a = BXYZ.ListOfBaseXyz()
    aStrXml = a.toStrXml()
    if verbose: print("test_210_0\n", aStrXml)
    self.assertEqual('<ListOfBaseXyz' in aStrXml, True)
    b = BXYZ.BaseFreeXyz()
    a.append(b)
    xml = a.toXml()
    #if verbose: print "test_210_1\n", ET.dump(xml)
    aStrXml = a.toStrXml()
    if verbose: print("test_210_2\n", aStrXml)
    self.assertEqual('<ListOfBaseXyz' in aStrXml, True)
    self.assertEqual('<BaseFreeXyz' in aStrXml, True)
    b.toto = BXYZ.BaseFreeXyz()
    b.tutu = BXYZ.ListOfBaseXyz()
    aStrXml = a.toStrXml()
    #if verbose: print "test_210_3\n", aStrXml
    """
    something like:
    <?xml version='1.0' encoding='UTF-8'?>
    <ListOfBaseFreeXyz size='1' typeClass='ListOfBaseXyz'>
      <BaseFreeXyz typeClass='BaseFreeXyz' index='0'>
        <toto typeClass='BaseFreeXyz'/>
        <tutu size='0' typeClass='ListOfBaseXyz'/>
      </BaseFreeXyz>
    </ListOfBaseFreeXyz>
    """
    self.assertEqual('<ListOfBaseXyz' in a.toStrXml(), True)
    self.assertEqual('<BaseFreeXyz' in a.toStrXml(), True)  

  def test_220(self):
    a = BXYZ.ListOfBaseXyz()
    b = BXYZ.BaseFreeXyz()
    c = BXYZ.BaseFreeXyz()
    a.append(b)
    a[0].toto = c
    aStrXml = a.toStrXml()
    if verbose: print("test_220_0\n", aStrXml)
    
    self.assertEqual(a.getNameAsAttribute() , None)
    self.assertEqual(b.getNameAsAttribute() , None) #is not an attribute of a, is element of list a
    self.assertEqual(a.getPyName() , "")
    self.assertEqual(b.getPyName() , "[0]")
    self.assertEqual(c.getPyName() , ".toto")
    self.assertEqual(c.getTreePyName() , "[0].toto")
    
    root = BXYZ.BaseFreeXyz()
    root.aa = a
    self.assertEqual(c.getTreePyName() , ".aa[0].toto")
    self.assertEqual(b.getTreePyName() , ".aa[0]")
    
    d = BXYZ.ListOfBaseXyz()
    e1 = BXYZ.BaseFreeXyz()
    e2 = BXYZ.BaseFreeXyz()
    d.append(e1)
    d.append(e2)
    self.assertEqual(e1.getNameAsAttribute() , None)
    self.assertEqual(e2.getNameAsAttribute() , None)
    self.assertEqual(e1.parentAsAttribute(), d)
    self.assertEqual(e2.parentAsAttribute(), d)
    self.assertEqual(d.getPyName() , "")
    self.assertEqual(e1.getPyName() , "[0]")
    self.assertEqual(e2.getPyName() , "[1]")
    self.assertEqual(d.getTreePyName() , "")
    self.assertEqual(e1.getTreePyName() , "[0]")
    self.assertEqual(e2.getTreePyName() , "[1]")
    
    root.bb = d
    self.assertEqual(d.getTreePyName() , ".bb")
    self.assertEqual(e1.getTreePyName() , ".bb[0]")
    self.assertEqual(e2.getTreePyName() , ".bb[1]")
    aStrXml = root.toStrXml()
    if verbose: print("test_220_1\n", aStrXml)
    
    """
    something like:
    <?xml version='1.0' encoding='UTF-8'?>
    <BaseFreeXyz typeClass='BaseFreeXyz'>
      <aa size='1' typeClass='ListOfBaseXyz'>
        <BaseFreeXyz typeClass='BaseFreeXyz' index='0'>
          <toto typeClass='BaseFreeXyz'/>
        </BaseFreeXyz>
      </aa>
      <bb size='2' typeClass='ListOfBaseXyz'>
        <BaseFreeXyz typeClass='BaseFreeXyz' index='0'/>
        <BaseFreeXyz typeClass='BaseFreeXyz' index='1'/>
      </bb>
    </BaseFreeXyz>
    """

  def test_300(self):
    # only on verbose human eye
    if not verbose: return
    s = """<?xml version='1.0' encoding='UTF-8'?>
<AROOTNAME typeClass='BaseFreeXyz'>
  <aa size='1' typeClass='ListOfBaseXyz'>
    <BaseFreeXyz typeClass='BaseFreeXyz' index='0'>
      <toto typeClass='BaseFreeXyz'/>
    </BaseFreeXyz>
  </aa>
  <bb size='2' typeClass='ListOfBaseXyz'>
    <BaseFreeXyz typeClass='BaseFreeXyz' index='0'/>
    <BaseFreeXyz typeClass='BaseFreeXyz' index='1'/>
  </bb>
</AROOTNAME>"""
    
    """bfree = BXYZ.BaseFreeXyz()
    #bxml = ET.fromstring(s)
    #bstr = UXYZ.prettyPrintET(bxml)
    #print "bstr",bstr
    bfree.fromXml(s)
    print "bfree from xml:\n",bfree
    print "bfree to xml:\n",bfree.toStrXml()"""
    bfree = UXYZ.fromXml(s)
    #bfreeStrXml = (bfree.toStrXml()).replace('"',"'").replace('utf-8',"UTF-8")
    bfreeStrXml = bfree.toStrXml()
    if verbose: print("!!!! compare initial and final:\n\ninitial:\n%s\n\nfinal:\n%s"%(s,bfreeStrXml))
    #self.failUnlessEqual(s.split("\n"), bfreeStrXml.split("\n"))

  def test_400(self):
    import xyzpy.intFloatListXyz as IFLX
    a = BXYZ.BaseFreeXyz()
    a.abasefree = BXYZ.BaseFreeXyz()
    a.int = IFLX.IntXyz(123)
    a.float = IFLX.FloatXyz(456.)
    a.str = IFLX.StrXyz("one Str")
    aStrXml = a.toStrXml()
    if verbose: print("\ntest_400_1\n", aStrXml)
    self.assertEqual( "<abasefree " in aStrXml, True)
    self.assertEqual( "<int " in aStrXml, True)
    self.assertEqual( "<float " in aStrXml, True)
    self.assertEqual( "<str " in aStrXml, True)
    
    self.assertEqual( "123</" in aStrXml, True)
    self.assertEqual( "one Str</" in aStrXml, True)
    self.assertEqual( "4" in aStrXml, True)

  """_XyzConstrainBase tests"""
  def test_510(self):
    a = BXYZ._XyzConstrainBase()
    self.assertEqual(a.isSet(), False)
    a.setIsSet(False) #useless for Free, useful for inherited
    self.assertEqual(a.isSet(), False)
    a.setIsSet(True)
    self.assertEqual(a.isSet(), True)
    
    self.assertEqual(a.getAttributes(), [])
    self.assertEqual(a.getCurrentAttributes(), [])
    self.assertRaises(Exception , a.__setattr__, "toto", 1 )
    
  def test_520(self):
    a = BXYZ._XyzConstrainBase()
    b = BXYZ._XyzConstrainBase()
    
    #it is a choice important do not change it
    #a == b only if ids are equals
    self.assertNotEqual(a, b) 
    self.assertEqual(a.equal(b), True)
    self.assertEqual(a.equal(None), False)
    self.assertRaises(Exception , a.equal,  1)
    
  def xtest_530(self):
    a = BXYZ._XyzConstrainBase()
    b = BXYZ._XyzConstrainBase()
    
    #a.toto = b #toto not expected
    self.assertRaises(Exception , a.__setattr__, "toto", b ) #a.toto = b
    
    self.assertEqual("toto" in list(a.__dict__.keys()), False)
    
    a._toto = b #accepted on name with previous '_'
    self.assertEqual("_toto" in list(a.__dict__.keys()), True)
    self.assertRaises(Exception , a.getCurrentAttributes) #not set
    self.assertRaises(Exception , a.getAttributes) #not set
    
    a.setIsSet(True)
    self.assertEqual(a.getCurrentAttributes(), []) #only on expected attributes
    self.assertEqual(a.getAttributes(), [])
        

  def test_540(self):
    
    class ExampleConstrainBaseXyz(BXYZ._XyzConstrainBase):
      _attributesList = [("toto", "StrXyz"), ("tutu", "IntXyz")]
      pass

    a = ExampleConstrainBaseXyz()
    b = ExampleConstrainBaseXyz()
    
    #a.tata = b #tata not expected
    self.assertRaises(Exception , a.__setattr__, "tata", b ) #a.tata = b
    #unexpected class ExampleConstrainBaseXyz for attribute toto, setIsCast False: no cast for the moment
    self.assertRaises(Exception , a.__setattr__, "toto", b ) #a.toto = b unexpected class ExampleConstrainBaseXyz for attribute toto
    
    a.setIsCast(True) 
    #setIsCast True: cast for now
    #b is castable to str with __str__ or __repr__
    a.toto = b #b is castable class to str, something as str "ExampleConstrainBaseXyz = [['_isSet', 'bool', False]]"
    
    a.toto = "hello"
    self.assertEqual("toto" in list(a.__dict__.keys()), True)
    self.assertEqual(a.toto, "hello")
    
    self.assertEqual("tutu" in list(a.__dict__.keys()), False)
    a.tutu = 1
    self.assertEqual("tutu" in list(a.__dict__.keys()), True)
    self.assertEqual(a.tutu, 1) #cast to IntXyz
    a.toto = 1 #cast to str
    self.assertEqual(a.toto, "1")
    a.toto = -01.2e-010
    self.assertEqual(a.toto, "-1.2e-10") #cast to str


if __name__ == '__main__':
  unittest.main()
  pass
