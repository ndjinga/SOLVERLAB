#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
utilities for and with ElementTree
InstanceToXml is a small browser to explore an instance of an objet python
recursively and obtain an elementTree (and a xml file) (for debugging)
"""

import sys
import xml.etree.ElementTree as ET
#from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem
import xyzpy.utilsXyz as UXYZ
import traceback 

#from copy import deepcopy
from PyQt5 import QtGui, QtWidgets

verbose = False
warningsDone = [] #store warnings for only one print

def printTrace():
   """use of traceback
   http://quentel.pierre.free.fr/python-trad/module-sys.html
   Avertissement: affecter la valeur de retour de trace à une variable locale dans une fonction
   provoquera une référence circulaire. 
   Ceci empêchera tout objet référencé par une variable locale 
   dans la même fonction ou dans la trace d'être nettoyé par le ramasse-miettes. 
   Puisque la plupart des fonctions n'ont pas besoin d'accéder à la trace, 
   la meilleure solution est d'utiliser quelque chose comme type, valeur = sys.exc_info()[:2] 
   pour extraire seulement le type d'exception et la valeur. 
   Si vous avez vraiment besoin de la trace, assurez-vous de la détruire après usage 
   (le meilleur moyen est d'utiliser une instruction try ... finally) 
   ou d'appeler exc_info() dans une fonction qui ne gère pas elle-même d'exception.
   """
   exctype, value, trace = sys.exc_info()
   traceback.print_exception(exctype, value, trace)
   
def printDir(obj,title=""):
    """a more detailed 'for i in dir(obj) : print i'"""
    check=[["Classes", "(class ", [], "class"],
           ["MethodWrappers", "(method-wrapper ", [], "MethodWrapper"],
           ["BoundMethods", "(bound method ", [], "BoundMethod"],
           ["BuiltInMethods", "(built-in method", [], "BuiltInMethod"],
           ["Objects", " object at ", [], "Object"],
           ["Instances", " instance at ", [], "Instance"],
           ["Others", "", [], "Other"]
          ]
    for dd in dir(obj):
      text=[dd, str(myEvalObj(obj, dd)).replace('<','(').replace('>',')')]
      for cc in check:
        if cc[1] in text[1]:
          cc[2].append(text)
          break
    print("\n---",title," ",obj.__class__)
    for cc in check:
      print(cc[0])
      for ii in cc[2]: print("   ",ii)
  
class BrowserPythonClass(object):
  """
  a small browser of python instance object (for debug)
  
  | example of use: 
  | explore methods etc. of a variable, variable level deep.
  | could save result in a file xml, or display it!
  | 
  | example: 
  |  >>> deepOfExplore = 4
  |  >>> res = BrowserPythonClass(aVariable, deepOfExplore)
  |  >>> res.toFileXml("myfile.xml")
  |  >>> res.display()
  """
  def __init__(self, obj, levelMax=2):
    #body=ET.Element("INSTANCE")
    body=ET.Element(obj.__class__.__name__)
    self.explore(body, obj, levelMax)
    self.body=body
  
  def strClass(self, obj):
    """if <class 'XFlux.XFluxOrg'> return XFluxOrg"""
    try: 
      res=str(obj.__class__).split('.')[-1]
      res=res.split("'")[0]
      return res
    except:
      print('WARNING: strClass: problem on %s.__class__' % str(obj))
      return str(type(obj))

  def strCompleteClass(self, app):
    """if <class 'XFlux.XFluxOrg'> return XFlux.XFluxOrg"""
    try: 
      res=str(app.__class__).split("'")[-2]
      return res
    except:
      print('WARNING: strCompleteClass: problem on %s.__class__' % str(obj))
      return str(type(obj))
  
  def strForXmlTag(self,  tag):
    """
    filter tag name as integer or float incorrect xml syntxx
    not well-formed (invalid token) for minidom.parseString
    """
    if issubclass(tag.__class__,  (int,  float)):
      return "_%s_" % str(tag)
    return str(tag)
    
  def toStrXml(self):
    """outputs exploration tree as xml string"""
    return UXYZ.prettyPrintET(self.body)
    
  def toFileXml(self, fileName, **kwargs):
    """outputs exploration tree as xml file"""
    strXml = self.toStrXml(**kwargs)
    file = open(fileName, "w")
    file.write(strXml)
    file.close()
    
  def display(self):
    """display with TreeXmlXyz xml object browser"""
    #app = OnceQApplication()
    self.treeWid = TreeXmlXyz()
    self.treeWid.setFromXml(self.body)
    self.treeWid.resize(500, 500)
    self.treeWid.show()
    return self.treeWid
  
  def _evalObj(self,  obj, string, default=""):
    """try eval(string) and value by default"""
    try:
      attr = getattr(obj, string)
      res = str(attr)
      #res = getattr(obj, string).__repr__()
    except:
      mess = 'WARNING: problem on _evalObj(%s.%s) : %s' % (type(obj), string, str(sys.exc_info()[1]))
      #print obj
      if mess not in warningsDone:
        warningsDone.append(mess)
        print(mess)
      res = 'problem on _evalObj(%s.%s)' % (type(obj), string)
    return res
    
  def _beginWith(self, aString, beginsWith ):
    lg = len(beginsWith)
    if aString[0:lg] == beginsWith:
      return True
    else:
      return False
      
  def explore(self, body, obj, levelMax=1, level=0):
    """"
    recursive exploration of a python object
    result is an Xml ElementTree
    """
    if level > levelMax: return
    levelcou = level+1
    check=[["Class", "(class ", [], "class"],
           ["MethodWrappers", "(method-wrapper ", [], "MethodWrapper"],
           ["BoundMethods", "(bound method ", [], "BoundMethod"],
           ["BuiltInMethods", "(built-in method", [], "BuiltInMethod"],
           #["Objects", " object at ", [], "Object"],
           #["Instances", " instance at ", [], "Instance"],
           ["Others", "", [], "Other"]
          ]
    
    thereIsClass = False
    skip = 0
    #test short exploration for elementary classes, if obj.__class__ exists
    try:
      classe = obj.__class__
      thereIsClass = True
    except:
      print("type %s has no attribute '__class__'" % type(obj))
      body.text = str(obj)
      body.attrib = {"type": str(type(obj))}
    
    if thereIsClass:
      if classe in [float, int,  str, None.__class__,  bool,  str]:
        body.text = str(obj)
        body.attrib = {"type": obj.__class__.__name__}
        return #shortest description for elementaries types
      
      skip = 1 #long description by default
      if classe in [list, tuple]:
        body.text = str(obj)
        body.attrib = {"type": obj.__class__.__name__}
        skip = 4 #short description by default
      
      if classe in [dict]:
        body.text = str(obj)
        body.attrib = {"type": obj.__class__.__name__}
        skip = 4 #short description by default
    
    #need long description...
    #[general type , _evalObj begins, [dd in dir() found], elementary text definition]
    check=[["Class", "(class ", [], "class"], #not to be in "Others"
           ["MethodWrappers", "(method-wrapper ", [], "MethodWrapper"],
           ["BoundMethods", "(bound method ", [], "BoundMethod"],
           ["BuiltInMethods", "(built-in method", [], "BuiltInMethod"],
           #["Objects", " object at ", [], "Object"],
           #["Instances", " instance at ", [], "Instance"],
           ["Others", "", [], "Other"]
          ]
          
    for dd in dir(obj):
      #avoid < and > in xml...
      text = [dd, str(self._evalObj(obj, dd)).replace('<','(').replace('>',')')]
      if verbose: print(text)
      for cc in check:
        if self._beginWith(text[1], cc[1]):
          cc[2].append(text)
          break
          
    for cc in check[skip:]: #skip class and standard method of elementary types(for economy)
      ii = ET.SubElement(body, cc[0])
      for idir, evalObj in cc[2]:
        jj = ET.SubElement(ii, idir)
        if cc[3] == "MethodWrapper":
          jj.text = getattr(obj, idir).__doc__
        elif cc[3] == "BoundMethod":
          jj.text = getattr(obj, idir).__doc__
        elif cc[3] == "BuiltInMethod":
          jj.text = getattr(obj, idir).__doc__
        elif cc[3] == "Other":
          try:
            attr = getattr(obj, idir)
          except:
            exctype, value, trace = sys.exc_info()
            jj.text = str(exctype)+str(value)
            continue
          #print idir, attr
          if idir == "__dict__":
            #sometimes cause problem jj.text = str(attr)
            jj.text = "keys = " + str(sorted(attr.keys()))
            #jj.attrib = {"keys": sorted(attr.keys())}
            continue
          if idir == "__doc__":
            jj.text = str(attr)
            continue
          if idir == "__module__":
            jj.text = str(attr)
            continue
          if idir == "__weakref__":
            jj.text = str(attr)
            continue
          if idir == "__class__":
            jj.text = str(attr)
            continue
          
          jj.text = str(attr)
          self.explore(jj, attr, levelMax, level=levelcou)
          continue
    
    if thereIsClass:
      if issubclass(classe,  dict):
        #jj.text = str(attr)
        #jj.attrib = {"keys": sorted(attr.keys())}
        nn = ET.SubElement(body, "dictItems")
        for key, value in list(obj.items()):
          kk = ET.SubElement(nn, self.strForXmlTag(key))
          kk.text = str(value)
          self.explore(kk, value, levelMax, level=levelcou)
        
      if issubclass(classe,  (list,  tuple)):
        #jj.text = str(attr)
        #jj.attrib = {"length": len(attr)}
        idx = 0
        nn = ET.SubElement(body, "ListItems")
        for value in obj:
          kk = ET.SubElement(nn, "_%i_" % idx)
          idx += 1
          kk.text = str(value)
          self.explore(kk, value, levelMax, level=levelcou)

