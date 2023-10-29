#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
utilities for using documentation sphinx in salome

example of use:

import pprint as PP
import utilSphinxSalome as USS

UT = USS.UtilSphinxSalome()
print(UT.getToday())
UT.setlocale()
print(UT.getToday())
inv = UT.fetch_inventory('file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/')
print(UT.get_c_locale_abbrev())
PP.pprint(inv.keys())
PP.pprint(inv)
inv = UT.fetch_inventory('file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/PACKAGESPY/doc/')
PP.pprint(inv.keys())
"""


import sys
import os
#import distutils.core
#import glob
#import shutil
#import datetime
#import subprocess
import fnmatch
import warnings
import pprint as PP
import time
import locale

from sphinx.ext.intersphinx import fetch_inventory

print("*** PYTHON version utilSphinxSalome.py", sys.version)

verbose = True

#######################################################################################
class UtilSphinxSalome(object):
  
  """utilities sphinx in salome, method of class used in files conf.py for sphinx."""
  
  def __init__(self):
    if verbose:
      print("UtilSphinxSalome.__init__")
    pass
    
  def isUnderCMake(self):
    """test from environment"""
    raise Exception("TODO isUnderCMake")
    return True
    
  def isUnderSat(self):
    """test from environment"""
    raise Exception("TODO isUnderSat")
    return True
  
  def isUnderSalome(self):
    """test from environment salome SalomeAppConfig set"""
    res = os.getenv("SalomeAppConfig")
    return (res != None)
    
  def getSalomeAppConfig(self):
    res = os.getenv("SalomeAppConfig")
    if res == None: return None
    return res.split(":")
  
  def getSalomePath(self):
    """get environment variable SALOMEPATH (existing under appli salome shell)"""
    """
    #example:
    SALOMEPATH=
      /export/home/salome/V7_5_1_CO6.4_64/INSTALL/KERNEL:
      /export/home/salome/V7_5_1_CO6.4_64/INSTALL/GUI:
      /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MATIX_PROFILE:
      /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MICROGEN:
      etc...
    """
    res = os.getenv("SALOMEPATH")
    if res == None: return None
    return res.split(":")
  
  def getRootDirModules(self,modules):
    res = []
    for mo in modules:
      var = os.getenv("%s_ROOT_DIR" % mo)
      print("%s_ROOT_DIR=%s" % (mo,var))
      if var != None: res.append(var)
    if len(res) == 0:
      return None
    else:
      return res
  
  def get_environs(self, search=""):
    env = os.environ #a dict
    res = dict((k,env[k]) for k in list(env.keys()) if search in k)
    if verbose:
      print("get_environs '%s'" % search)
      pprint.pprint(res)
    return res
    
  def print_environment(self):
    print("\nenvironment:\n%s" % PP.pformat(sorted(os.environ.keys())))
  
  def recursive_glob(self, treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
      goodfiles = fnmatch.filter(files, pattern)
      results.extend(os.path.join(base, f) for f in goodfiles)
    return results
    
  def find_file(self, rootPath, name):
    """find first file named in root recursively"""
    res = self.recursive_glob(rootPath, name)
    if len(res) == 0: return None
    return res[0]
   
  def getInstalledDoc(self, rootPath):
    """find first entry point index.html of an installed sphinx doc, from root file path"""
    return self.find_file(rootPath, "index.html")
  
  def get_intersphinx_data(self, modules):
    """
    get paths to modules docs salome from environment
    """
    """
    look for ? TODO ?
    SRC_ROOT=/export/home/wambeke/MATIX_V21_S751_CO6.4_64/SOURCES

    DOCUTILS_ROOT_DIR
    CRESCENDO_ROOT_DIR
    CRESCENDO_SRC_DIR
    FFTWCODE_ROOT_DIR
    FFTWCODE_SRC_DIR
    INSTALL_ROOT=/export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL
    MATIX_PROFILE_ROOT_DIR=/export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MATIX_PROFILE
    MATIX_PROFILE_SRC_DIR=/export/home/wambeke/MATIX_V21_S751_CO6.4_64/SOURCES/MATIX_PROFILE
    PRODUCT_ROOT_DIR=/export/home/wambeke/sphinxpy
    SALOMEPATH=
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/KERNEL
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/GUI
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MATIX_PROFILE
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MICROGEN
  etc...
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/PACKAGESPY
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/CODES/FFTWCODE
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/CODES/MFRONTCODE
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/CODES/AMITEXCODE
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/LIBBATCH
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/GEOM
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/SMESH
  ...
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/JOBMANAGER
    SALOME_APPLI_ROOT=/export/home/wambeke/MATIX_V21_S751_CO6.4_64/APPLI
    SALOME_MODULES=MATIX_PROFILE,MICROGEN,CRESCENDO,NUMODIS,AMITEX,CMDC,EKINOX,DART,IRADINAGUI,PACKAGESPY,FFTWCODE,MFRONTCODE,AMITEXCODE,LIBBATCH,KERNEL,GUI,GEOM,SMESH,MED,YACS,NETGENPLUGIN,GHS3DPLUGIN,XDATA,YACSGEN,PARAVIS,HEXABLOCK,HEXABLOCKPLUGIN,JOBMANAGER
    SalomeAppConfig=
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MATIX_PROFILE/share/salome/resources/matix_profile
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/GUI/share/salome/resources/gui
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MATIX_PROFILE/share/salome/resources/matix_profile
  /export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/MICROGEN/share/salome/resources/microgen
  etc...
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/HEXABLOCKPLUGIN/share/salome/resources/hexablockplugin
  /export/home/salome/V7_5_1_CO6.4_64/INSTALL/GHS3DPLUGIN/share/salome/resources/ghs3dplugin
    """
    rst_prolog = ""
    intersphinx_mapping = {}
    links_rst = """\
.. sphinxpy all documentation rst file, created by
   UtilSphinxSalome.py method get_intersphinx_data()
   !!!This file is overriden in make html!!!


.. _entryalldocumentations:


Toutes documentations
======================

Cette page vous donne accès aux documentations des modules présents.

%s

"""

    toctree = ""
    #envs = self.get_environs("INSTALL")
    modules_path = {}
    
    #try from SalomeAppConfig
    #modules_install = self.getSalomeAppConfig()
    #try from SALOMEPATH
    modules_install = self.getSalomePath()
    if modules_install == None:
      modules_install = self.getRootDirModules(modules)
    
    print("modules_install %s" % modules_install)
    if modules_install == None:
      print("\nproblem SalomePath not found")
      self.print_environment()
      modules_install = []
      
    for mo in modules:
      for pa in modules_install:
        if mo in pa:
          doc = self.getInstalledDoc(pa)
          modules_path[mo] = [pa, doc]
          break

    rst = ""
    #for mo, (pa, doc) in modules_path.items():
    for mo in modules:
      try:
        pa, doc = modules_path[mo]
      except:
        continue
      aFile = """file://%s""" % doc
      rst_prolog += """.. _%s: %s\n""" % (mo, aFile)
      intersphinx_mapping[mo] = aFile
      rst += " * La documentation du module %s_\n" % mo
      
    links_rst = links_rst % rst
    
    if verbose: 
      print("\nget_intersphinx_data:\n%s" % PP.pformat(intersphinx_mapping))
      #pprint.pprint(modules_path)
      print("\ntoutesdocumentations.rst:\n\n%s\n" % links_rst)
    
    return rst_prolog, intersphinx_mapping, links_rst
    

  def get_c_locale_abbrev(self):
    """for example, get english time even if french locale set"""
    lc = locale.setlocale(locale.LC_TIME) #initial locale
    try:
      locale.setlocale(locale.LC_TIME, "C") #english locale
      res = time.strftime("%a-%d-%b") #'Thu-05-Nov'
    except:
      res = time.strftime("%a-%d-%b") #'jeu.-05-nov.' in french
    #set initial locale
    locale.setlocale(locale.LC_TIME, lc) 
    return res

  def setlocale(self, value='fr_FR.utf8'):
    """suppose that we want to be french by default"""
    locale.setlocale(locale.LC_ALL, 'fr_FR.utf8')
    
  def getToday(self):
    """needs setlocale fr done before to have french date"""
    #res = unicode(time.strftime('%A %d %B %Y')) #mercredi 21 janvier 2015 or 05 November 2015
    #fix french accent problem if "février"res = str(time.strftime('%d %B %Y').decode('utf-8')) #fix french accent problem if "février"
    # obsolete python2 res = u"%s" % time.strftime('%d %B %Y').decode('utf-8')

    # python3
    b = bytes(time.strftime('%d %B %Y'), 'utf-8')
    res = str(b, 'utf-8')
    return res

  def fetch_inventory(self, uri):
    """
    Read inventory into a dictionnary to explore sphinx file objects.inv.
    uri could be local file path (warning for file needs three '/'):
    
    - 'http://docs.python.org/2.7/'
    - 'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/'
    - etc...
    """

    #http://stackoverflow.com/questions/30939867/how-to-properly-write-coss-references-to-external-documentation-with-intersphinx 

    if verbose: print("fetch_inventory %s" % uri)
    inv = fetch_inventory(warnings, uri, uri + 'objects.inv')
    if verbose: 
      print("fetch_inventory keys:\n%s" % PP.pformat(list(inv.keys())))
    return inv

    """
uri = 'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/'

# Read inventory into a dictionnary
inv = fetch_inventory(warnings, uri, uri + 'objects.inv')
pprint.pprint(inv)
{u'std:label': {u'genindex': (u'numodis',
                              u'2.1',
                              u'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/genindex.html#',
                              u'Index'),
                u'modindex': (u'numodis',
                              u'2.1',
                              u'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/py-modindex.html#',
                              u'Module Index'),
                u'search': (u'numodis',
                            u'2.1',
                            u'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/NUMODIS/share/doc/salome/gui/numodis/search.html#',
                            u'Search Page')}}

uri = 'file:///export/home/wambeke/MATIX_V21_S751_CO6.4_64/INSTALL/PACKAGESPY/doc/'

# Read inventory into a dictionnary
inv = fetch_inventory(warnings, uri, uri + 'objects.inv')
pprint.pprint(inv.keys())
[u'std:label',
 u'py:class',
 u'py:function',
 u'py:method',
 u'py:module',
 u'py:data',
 u'py:attribute']
"""

#######################################################################################
if __name__ == '__main__':
  if verbose:
    print(__doc__)
  pass
  





