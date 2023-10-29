#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import unittest
from xyzpy.baseXyz import BaseFreeXyz, ListOfBaseXyz
import xyzpy.intFloatListXyz as IFLX

def root_path():
  return os.path.abspath(os.sep)

verbose=False

###############################################################
class essaiSetattr(object):
  # elementary example
  def __init__(self):
    self.__dict__['_paramOrder'] = []

  def __setattr__(self, name, value):
    if verbose: print("__setattr__", name, value)
    if name in self._paramOrder:
      self.__delattr__(name)
    self._paramOrder.append(name)
    object.__setattr__(self, name, value)

  def __delattr__(self, name):
    if verbose: print("__delattr__", name, value)
    self._paramOrder.remove(name)
    object.__delattr__(self, name)

  def __iter__(self):
    for name in self._paramOrder:
      yield getattr(self, name)

  def items(self):
    for name in self._paramOrder:
      yield name, getattr(self, name)


###############################################################
class essaiSetattrBase(object):
  # more specialized example
  def __init__(self):
    self._attributes = [ #list, not a dict because sequential order list is used in files cnf
      ("bb", IFLX.FloatPosXyz),
      ("cc", int),
      ("dd", float),
      ("aa", IFLX.Int01Xyz),
    ]
    self._attributesDict=dict((key, value) for (key, value) in self._attributes)

  def __setattr__(self, name, value):
    if verbose: print("__setattr__", name, value)
    if name[0] == "_":
      return object.__setattr__(self, name, value)
    if name in list(self._attributesDict.keys()):
      if name in list(self.__dict__.keys()):
        self.__delattr__(name)
      obj = self._attributesDict[name](value)
      return object.__setattr__(self, name, obj)
    else:
      raise Exception(self.__class__.__name__ + ": attribute '" + name + "' is unexpected")

  """def __delattr__(self, name):
    if verbose: print "__delattr__", name, value
    self._paramOrder.remove(name)
    object.__delattr__(self, name)"""

  def __iter__(self):
    for name, _ in self._attributes:
      if name in list(self.__dict__.keys()):
        yield getattr(self, name)
      else:
        raise Exception(self.__class__.__name__ + ": attribute '" + name + "' is undefined")

  def items(self):
    for name, _ in self._attributes:
      if name in list(self.__dict__.keys()):
        yield name, getattr(self, name)
      else:
        raise Exception(self.__class__.__name__ + ": attribute '" + name + "' is undefined")


###############################################################
class essaiSetattrBaseParentAsAttribute(object):
  # more specialized example, with parentAsAttribute() function to get parent as in tree for attributes
  def __init__(self):
    self._attributes = [ #list, not a dict because sequential order list is used in files cnf
      ("bb", IFLX.FloatPosXyz),
      ("cc", IFLX.IntXyz),
      ("dd", IFLX.FloatXyz),
      ("aa", IFLX.Int01Xyz),
      ("aTree", essaiSetattrBaseParentAsAttribute),
    ]
    self._attributesDict=dict((key, value) for (key, value) in self._attributes)
    #only one modification allowed, only one parent in tree, no loop.
    #set in instance value in __setattr__() of parent
    self._parentAsAttribute = None

  def __setattr__(self, name, value):
    if verbose: print("__setattr__", name, value)

    if name[0] == "_": return object.__setattr__(self, name, value) #a classical attribute

    if name not in list(self._attributesDict.keys()): raise Exception(self.__class__.__name__ + ".__setattr__: attribute '" + name + "' is unexpected")

    obj = value #default
    if hasattr(value, "_parentAsAttribute"):
      if value._parentAsAttribute != None:
        raise Exception(self.__class__.__name__ + ".__setattr__: parentAsAttribute of value is set yet")

    else: #value without _parentAsAttribute: integer for non mutable Int (for example)
      if True: #try:
        obj = self._attributesDict[name](value) #obj=Int(value) (for example)
        if verbose: print('create obj as attribute:',name,self._attributesDict[name],obj)
      else: #except:
        raise Exception(self.__class__.__name__ + 
          ".__setattr__: value incompatible with conversion to " + 
          self._attributesDict[name].__name__ + 
          ": " + str(value))

    if name in list(self.__dict__.keys()): #attribute existing yet
        self.__delattr__(name)

    object.__setattr__(self, name, obj)
    obj._parentAsAttribute = self
    return

    """ 
    #thanks to introspection, for obj class without _parentAsAttribute mechanism
    #may be for instances: obj._parentAsAttribute = self
    #could be risky: pass
    raise Exception(obj.__class__.__name__ + ".__setattr__: value without _parentAsAttribute attribute")
    """
  def __iter__(self):
    for name, _ in self._attributes:
      if name in list(self.__dict__.keys()):
        yield getattr(self, name)
      else:
        raise Exception(self.__class__.__name__ + ": attribute '" + name + "' is undefined")

  def items(self):
    for name, _ in self._attributes:
      if name in list(self.__dict__.keys()):
        yield name, getattr(self, name)
      else:
        raise Exception(self.__class__.__name__ + ": attribute '" + name + "' is undefined")

  def parentAsAttribute(self):
    return self._parentAsAttribute


###############################################################
class TestCase(unittest.TestCase):
  def iterAll(self, a):
    return [i for i in a]
  
  def iteritemAll(self, a):
    return [(name, val) for name, val in list(a.items())]

  def setattribute(self, a, name, value):
    a.__setattr__(name, value)
  
  def test_010(self):
    # elementary example
    a=essaiSetattr()
    self.assertEqual(a.__dict__, {'_paramOrder': []})
    self.assertEqual(a._paramOrder, [])
    a.toto = 1
    self.assertEqual(a.toto, 1)
    self.assertEqual(sorted(a.__dict__.keys()), ['_paramOrder', 'toto'])
    a.tutu = 2
    self.assertEqual(a.tutu, 2)
    self.assertEqual(sorted(a.__dict__.keys()), ['_paramOrder', 'toto', 'tutu'])
    self.assertEqual(a._paramOrder, ['toto', 'tutu'])
    a.tata = 3
    self.assertEqual([i for i in a], [1, 2, 3])
    self.assertEqual([(name, val) for name, val in list(a.items())], [('toto', 1), ('tutu', 2), ('tata', 3)])
    del(a.tutu)
    self.assertEqual([(name, val) for name, val in list(a.items())], [('toto', 1), ('tata', 3)])
    self.assertEqual(a._paramOrder, ['toto', 'tata'])
    a.tutu = 4
    self.assertEqual([(name, val) for name, val in list(a.items())], [('toto', 1), ('tata', 3), ('tutu', 4)])
    a.tata = 5
    self.assertEqual([(name, val) for name, val in list(a.items())], [('toto', 1), ('tutu', 4), ('tata', 5)])
    self.assertEqual(sorted(a.__dict__.keys()), ['_paramOrder', 'tata', 'toto', 'tutu'])
    
  def test_020(self):
    # more specialized example
    a=essaiSetattrBase()
    self.assertEqual(sorted(a.__dict__.keys()), ['_attributes', '_attributesDict'])
    self.assertRaises(Exception , self.iterAll, a)
    self.assertRaises(Exception , self.iteritemAll, a)
    
    self.assertRaises(Exception , self.setattribute, a, "aa", 4)
    self.assertRaises(Exception , self.setattribute, a, "bb", -1)
    self.assertRaises(Exception , self.setattribute, a, "toto", 1)
    
    a.aa=1
    a.bb=2
    a.cc=3
    a.dd=4
    
    self.assertEqual(sorted(a.__dict__.keys()), ['_attributes', '_attributesDict', 'aa', 'bb', 'cc', 'dd'])
    self.assertEqual([val for val in a], [2.0, 3, 4.0, 1])
    self.assertEqual([(name, val) for name, val in list(a.items())], [('bb', IFLX.FloatPosXyz(2.0)), ('cc', 3), ('dd', 4.0), ('aa', IFLX.Int01Xyz(1))])
    
    a.aa=0
    self.assertEqual([val for val in a], [2.0, 3, 4.0, 0])
    
  def test_030(self):
    # more specialized example, with parentAsAttribute()
    a = IFLX.IntXyz(123)
    self.assertEqual('_parentAsAttribute' in list(a.__dict__.keys()), True)
    self.assertEqual(a.parentAsAttribute(), None)
    
    a = essaiSetattrBaseParentAsAttribute()
    self.assertEqual(sorted(a.__dict__.keys()), ['_attributes', '_attributesDict', '_parentAsAttribute'])
    self.assertRaises(Exception , self.iterAll, a)
    self.assertRaises(Exception , self.iteritemAll, a)
    
    self.assertRaises(Exception , self.setattribute, a, "aa", 4)
    self.assertRaises(Exception , self.setattribute, a, "bb", -1)
    self.assertRaises(Exception , self.setattribute, a, "toto", 1)
    
    a.aa=1
    a.bb=2
    a.cc=3
    a.dd=4
    a.aa=0
    b = essaiSetattrBaseParentAsAttribute()
    b.aa=0
    b.bb=22
    b.cc=33
    b.dd=44
    self.assertRaises(Exception , a.__setattr__, "aTree", 1234)
    
    a.aTree = b
    self.assertRaises(Exception , a.__setattr__, "aTree", b) #only once setattr of b
    
    self.assertEqual(
      sorted(a.__dict__.keys()), 
      ['_attributes', '_attributesDict', '_parentAsAttribute', 'aTree', 'aa', 'bb', 'cc', 'dd'])
    self.assertEqual(a.aa, 0)
    self.assertEqual(a.bb, 2)
    self.assertEqual(a.cc, 3)
    self.assertEqual(a.dd, 4)
    self.assertEqual(a.aTree.aa, 0)
    self.assertEqual(a.aTree.bb, 22)
    self.assertEqual(a.aTree.cc, 33)
    self.assertEqual(a.aTree.dd, 44)

    self.assertEqual(
      [(name, val) for name, val in list(a.items())], 
      [('bb', IFLX.FloatPosXyz(2.0)),
       ('cc', IFLX.IntXyz(3)),
       ('dd', IFLX.FloatXyz(4.0)),
       ('aa', IFLX.Int01Xyz(0)),
       ('aTree', b)] )
    for name, val in list(a.items()):
      self.assertEqual( val.parentAsAttribute(), a )

if __name__ == '__main__':
  unittest.main()
  pass
