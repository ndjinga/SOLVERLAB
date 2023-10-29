#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""\
tests for intFloatListXyz
"""

import os
import unittest
import xyzpy.intFloatListXyz as IFLX

verbose = False

class testInt01(int):
   def __new__(cls, value):
     """
     | see http://stackoverflow.com/questions/2673651/inheritance-from-str-or-int
     | see https://docs.python.org/2/reference/datamodel.html
     """
     if value in [0,1]:
       return int.__new__(cls, value)
     else:
       raise Exception("Int01: is not in [0,1]")
   
   def sameType(self, value):
     return self.__class__(value)

   def __repr__(self):
     #__repr__ goal is to be unambiguous
     cl = self.__class__
     return cl.__name__ + '(' + int.__repr__(self) + ')'

   def __str__(self):
     return int.__repr__(self)

   def sameType(self, value):
     return self.__class__(value)

class testIntRange(int):
   def __new__(cls, value, minMax):
     #print "new minMax",minMax
     if value >= minMax[0] and  value <= minMax[1]:
       obj = int.__new__(cls, value)
       obj.minMax = minMax
       return obj
     else:
       raise Exception("IntRange: " + str(value)+" is not in "+ str(minMax))
   
   def sameType(self, value):
     #print "sametype minMax", self.minMax
     return self.__class__(value, self.minMax)

   def __repr__(self):
     #__repr__ goal is to be unambiguous
     cl = self.__class__
     return cl.__name__ + '(' + int.__repr__(self) + ', ' + str(self.minMax) + ')'

   def sameType(self, value):
     return self.__class__(value, self.minMax)

class TestCase(unittest.TestCase):

  def test_010(self):
    a=0
    c=a
    idc=id(c)
    b=1
    self.assertNotEqual(a, b)
    self.assertNotEqual(id(a), id(b))
    self.assertEqual(id(a), id(c))
    b=0
    self.assertEqual(a, b)
    self.assertEqual(id(a), id(b)) #!!!(A) is here see below
    b=3
    self.assertNotEqual(a, b)
    self.assertNotEqual(id(a), id(b))
    a=b
    self.assertEqual(a, b)
    self.assertEqual(id(a), id(b))
    self.assertEqual(c, 0)
    self.assertEqual(id(c), idc)
    a=4
    self.assertNotEqual(a, b)
    self.assertNotEqual(id(a), id(b))
    self.assertEqual(b, 3)

  def test_020(self):
    a=testInt01(1)
    self.assertEqual(a, 1)
    self.assertEqual(str(a), "1")
    self.assertEqual(a.__repr__(), "testInt01(1)")
    a=testInt01(0)
    #for i in dir(a): print "dir a",i
    self.assertRaises(Exception , testInt01, -10)
    self.assertRaises(Exception , testInt01, 2)

    c=a
    b=testInt01(1)
    self.assertNotEqual(a, b)
    self.assertNotEqual(id(a), id(b))
    self.assertEqual(id(a), id(c))
    b=testInt01(0)
    self.assertEqual(a, b)
    self.assertNotEqual(id(a), id(b)) #!!!diff from (A)
    self.assertEqual(b+3, 3) #!!!diff from (A)
    self.assertEqual(isinstance(b,testInt01), True)
    self.assertEqual(isinstance(b+3, int), True)

    d=b.sameType(1)
    self.assertEqual(isinstance(d,testInt01), True)
    self.assertNotEqual(b, d)
    self.assertEqual(1, d)

  def test_030(self):
    minMax = [10,20]
    #a=IntRange(1, minMax=minMax)
    self.assertRaises(Exception , testInt01, 9, minMax)
    self.assertRaises(Exception , testInt01, 21, minMax)
    self.assertEqual(testIntRange(10, minMax), 10)
    self.assertEqual(testIntRange(11, minMax), 11)
    self.assertEqual(testIntRange(20, minMax), 20)
    self.assertEqual(testIntRange(20, minMax).__repr__(), "testIntRange(20, [10, 20])")

    a=testIntRange(10, [9,10])
    b=testIntRange(11, [11,12])
    self.assertRaises(Exception , a.sameType, 13)
    c=a.sameType(9)
    self.assertEqual(a, 10)
    self.assertEqual(c, 9)

    d=a.__class__(9, [9,90]) #cht range
    self.assertEqual(d, 9)
    self.assertEqual(c, d)
    self.assertEqual(c.__repr__(), "testIntRange(9, [9, 10])")
    self.assertEqual(d.__repr__(), "testIntRange(9, [9, 90])")

  def test_040(self):
  
    class TTEST(IFLX.StrInListXyz):
      _allowedList = ['hello','bye']
      pass
      
    a = TTEST('hello')
    self.assertEqual(a, "hello")
    self.assertRaises(Exception , TTEST, 'bughello')
    b = a.sameType('bye')
    self.assertEqual(b, "bye")
    self.assertRaises(Exception , a.sameType, 'bughello')
    
  def test_100(self):
    #for i in dir(IFLX): print i
    for key in list(IFLX.__dict__.keys()):
      if key[0] != "_":
        if key[-3:] != "Xyz": continue #not an immutable
        aClassType = IFLX.__dict__[key]
        if not issubclass(aClassType, IFLX._XyzImmBase): 
          if verbose: print("test intFloatList."+key+" not tested")
          continue
        if verbose: print("test intFloatList."+key)
        ii = 1 #default
        if "Release" in key: continue #test plus loin test_110
        if "Str" in key: continue #test plus loin test_115
        if "Bool" in key: continue #test plus loin test_125
        if "None" in key: continue #NoneXyz(None)
        if "Range" in key:
          a = aClassType(1)
        elif "InList" in key:
          a = aClassType(1)
        elif "Neg" in key:
          ii = -1
          a = aClassType(ii)
          self.assertRaises(Exception , aClassType, -ii)
        elif "Sup" in key:
          ii = 10
          a = aClassType(ii)
          self.assertRaises(Exception , aClassType, -ii)
        elif "Int" in key:
          ii = 1
          if "Int0Xyz" in key: ii = 0
          if "Int3Xyz" in key: ii = 3
          a = aClassType(ii)
          self.assertRaises(Exception , aClassType, "bug")
        elif "Float" in key:
          a = aClassType(ii)
          self.assertRaises(Exception , aClassType, "bug")
        else:
          try:
            a = aClassType(ii)
          except:
            raise Exception("Error: in %s unexpected"%key)
        #if verbose: print "not certain", a.__repr__()+" -> " + str(a)
        self.assertEqual(key in a.__repr__(), True)
        self.assertEqual(str(ii) in str(a), True)

  def test_110(self):
    #for i in dir(intFloatList): print i
    a = IFLX.ReleaseXyz()
    self.assertEqual("0.0.0" in a.__repr__(), True)
    self.assertEqual("0.0.0", str(a))
    self.assertRaises(Exception , IFLX.ReleaseXyz, "123.456")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "1.3.4.5")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "toto")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "-1")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "100")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "1.2.")
    self.assertRaises(Exception , IFLX.ReleaseXyz, "1..")

  def test_115(self):
    #for i in dir(intFloatListXyz): print i
    self.assertRaises(Exception , IFLX.StrModelXyz, "FeCu_P")
    a = IFLX.StrModelXyz("FeCu_p")
    #print a.__repr__()
    self.assertEqual(str(a), 'FeCu_p' )
    self.assertEqual(a.__repr__(), "StrModelXyz('FeCu_p', ['FeCu', 'FeCu_p', 'FeHe'])" )

  def test_120(self):
    self.assertRaises(Exception , IFLX.BoolWithIndeterminatedXyz, "oops")
    self.assertRaises(Exception , IFLX.BoolWithIndeterminatedXyz, "1")
    self.assertRaises(Exception , IFLX.BoolWithIndeterminatedXyz, "0")
    a = IFLX.BoolWithIndeterminatedXyz(True)
    #print a.__repr__()
    self.assertEqual(a.__repr__(), 'BoolWithIndeterminatedXyz(True)')
    self.assertEqual(a.getValue(), True)
    self.assertEqual(a, "True")
    self.assertEqual(bool(a), True)
    a = IFLX.BoolWithIndeterminatedXyz(False)
    self.assertEqual(a.getValue(), False)
    self.assertEqual(a, "False")
    self.assertEqual(a.__repr__(), 'BoolWithIndeterminatedXyz(False)')
    self.assertEqual(bool(a), False)
    a = IFLX.BoolWithIndeterminatedXyz()
    self.assertEqual(a.getValue(), None)
    self.assertEqual(a, "Indeterminated")
    self.assertEqual(a.__repr__(), 'BoolWithIndeterminatedXyz(None)')
    self.assertRaises(Exception , bool, a)
    
  def test_125(self):
    self.assertRaises(Exception , IFLX.BoolXyz, "oops")
    self.assertRaises(Exception , IFLX.BoolXyz, "1")
    self.assertRaises(Exception , IFLX.BoolXyz, "0")
    self.assertRaises(Exception , IFLX.BoolXyz, " ")
    self.assertRaises(Exception , IFLX.BoolXyz, "")
    
    a = IFLX.BoolXyz(True)
    self.assertEqual(a.getValue(), True)
    self.assertEqual(a, "True") #a is a str, not a bool (not inheritable)
    self.assertEqual(a.__repr__(), 'BoolXyz(True)')
    self.assertEqual(bool(a), True)
   
    a = IFLX.BoolXyz(False)
    self.assertEqual(a.getValue(), False)
    self.assertEqual(a, "False")
    self.assertEqual(a.__repr__(), 'BoolXyz(False)')
    self.assertEqual(bool(a), False)
    
    a = IFLX.BoolXyz()
    self.assertEqual(IFLX.BoolXyz(), "False") #is a str, not a bool not inheritable
    self.assertEqual(a.getValue(), False)
    self.assertEqual(a, "False") #a is a str, not a bool (not inheritable)
    self.assertEqual(a.__repr__(), 'BoolXyz(False)')
    self.assertEqual(bool(a), False)

  def test_130(self):
    self.assertEqual(IFLX.NoneXyz("oops"), "None")
    self.assertEqual(IFLX.NoneXyz(123), "None")
    
    a = IFLX.NoneXyz()
    self.assertEqual(a.getValue(), None)
    self.assertEqual(a, "None") #a is a str, not a None (not inheritable)
    self.assertEqual(a.__repr__(), 'NoneXyz(None)')
    self.assertEqual(bool(a), True) #it is a non empty str



if __name__ == '__main__':
  unittest.main()


