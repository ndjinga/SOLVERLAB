#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
utilities to WORKDIR set env and create directories
"""

import os
#import sys
import platform
import shutil
from time import sleep
from datetime import datetime
import subprocess
import pprint as PP #pretty print

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

_initializeWorkdirMatixDone = False
_workdirsMatix = ["DART", "CRESCENDO", "EKINOX", "MICROGEN", "NUMODIS", "AMITEX", "CMDC", "TESTS", "MATPLOTLIB", "IRADINAGUI"]

verbose = True


def getUserName():
  res = os.getenv('USERNAME')
  if res == None:
    res = os.getenv('USER')
  if res is None:
    raise Exception("can't get user name in env var 'USER' or 'USERNAME'")
  return res

def getHomeDir():
  res = os.getenv('HOME')
  if res == None:
    res = os.getenv('USERPROFILE')
  if res is None:
    raise Exception("can't get home dir in env var 'HOME' or 'USERPROFILE'")
  return res

def osCommand(cmd):
  if verbose: print("osCommand: '%s'" % cmd)
  process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  stdout, stderr = process.communicate()
  return (stdout, stderr)

def rmTree(aDir):
  #http://stackoverflow.com/questions/1557351/python-delete-non-empty-dir
  if aDir in ["", "/"]: raise Exception("rmTree : remove tree on '/' forbidden")
  shutil.rmtree(aDir, ignore_errors=True) #problemes si file system lent ???
  if os.path.exists(aDir): raise Exception("rmTree : remove tree NOT done: %s" % aDir)
  return
  '''
  #or ...
  if aDir in ["", "/"]: raise Exception("rmTree : remove tree on '/' forbidden")
  stdout, stderr = osCommand("rm -rf " + aDir)
  return stderr
  '''

def makeDir(nameDir):
  """
  if nameDir exists, make only one save copy,
  rename it as .bak, removing previous existing .bak
  """
  if os.path.exists(nameDir):
    dirbak=nameDir+".bak"
    if os.path.isdir(dirbak): 
      shutil.rmtree(dirbak)
    os.rename(nameDir,dirbak)
    os.listdir(dirbak) #sert seulement a mettre a jour le systeme de fichier sur certaines machines
  os.makedirs(nameDir)


def getWorkdirDefault(moduleName):
  """create working directory for moduleName if not existing"""
  try:
    workDir = os.getenv("DEFAULT_WORKDIR")
  except:
    workDir = None

  if workDir is None:
    logger.warning("environment var DEFAULT_WORKDIR not known, get ${HOME}/WORKDIR")
    workDir = os.path.join(os.getenv("HOME"), "WORKDIR")
    os.environ["DEFAULT_WORKDIR"] = workDir

  try:
    workDirModule = os.path.join(workDir, moduleName)
  except:
    logger.warning("DEFAULT_WORKDIR problem on moduleName '%s'" % moduleName)
    workDirModule = os.path.join(workDir, 'DEFAULT')

  if not os.path.exists(workDirModule):
    os.makedirs(workDirModule)
    logger.warning('create workdir ' + workDirModule)
  
  return workDirModule


def getWorkdirMatix(moduleName):
  if _initializeWorkdirMatixDone == False:
    initializeWorkdirMatix()
  try:
    return os.getenv(moduleName + "_WORKDIR")
  except:
    logger.error("%s_WORKDIR not known" % moduleName)
    res = getWorkdirDefault(moduleName)
    return res


def getWorkdir(moduleName):
  logger.warning("getWorkdir obsolescent method for MATIX '%s', use getWorkdirMatix()" % moduleName)
  return getWorkdirMatix(moduleName)


def obsolete_copyDir(srcDir, destDir):
  """
  if destDir exists, make only one save copy,
  rename it as .bak, removing previous existing .bak
  """
  if os.path.exists(destDir):
    dirbak=destDir + ".bak"
    if os.path.exists(dirbak):
      logger.debug("rmtree %s isdir=%s isfile=%s exists=%s" % \
                   (dirbak, os.path.isdir(dirbak), os.path.isfile(dirbak), os.path.exists(dirbak)))
      logger.debug("bak:\n%s" % PP.pformat(osCommand("ls -al " + dirbak)))
      #http://stackoverflow.com/questions/1557351/python-delete-non-empty-dir
      #shutil.rmtree(dirbak, ignore_errors=True) problemes si file system lent??
      rmTree(dirbak)
    os.rename(destDir, dirbak)
    os.listdir(dirbak) #sert seulement a mettre a jour le systeme de fichier sur certaines machines

  # shutil.copytree: The destination directory, named by dst, must not already exist
  # print "copyDir(srcDir, destDir)\n  %s\n  %s" % (srcDir, destDir), os.path.exists(destDir)
  rmTree(os.path.realpath(destDir))
  shutil.copytree(srcDir, destDir)

def copyDir(srcDir, destDir):
  """
  if destDir exists, delete it
  """
  realDestDir = os.path.realpath(destDir)
  # shutil.copytree: The destination directory, named by dst, must not already exist
  # print "copyDir(srcDir, realDestDir)\n  %s\n  %s" % (srcDir, realDestDir), os.path.exists(realDestDir)

  if os.path.exists(realDestDir):
    if os.path.isfile(realDestDir):
      logger.critical("copy directory to an existing file (not an existing directory):\n%s" % realDestDir)
    else:
      logger.warning("copy directory overriding existing directory:\n%s" % realDestDir)
    rmTree(realDestDir)
  shutil.copytree(srcDir, realDestDir)

def copyFile(src, dst, verb=False):
  if verb: print("create file: '%s'" % dst)
  shutil.copyfile(src, dst)
  return dst

def chmodarwx(nameFile):
  """chmod a+rwx"""
  import stat
  st = os.stat(nameFile)
  os.chmod(nameFile, st.st_mode | stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
  return
    
def initializeWorkdirMatix():
  """
  initialize user working directory $MATIX_WORKDIR if not existing
  """
  workDir = os.getenv("MATIX_WORKDIR")
  if workDir == None: #supposed not initialized at all
    os.environ["MATIX_WORKDIR"] = os.path.join(os.getenv("HOME"), 'MATIXWORKDIR')
    os.environ["DEFAULT_WORKDIR"] = os.path.join(os.getenv("HOME"), 'MATIXWORKDIR')
    workDirIni = os.getenv("MATIX_WORKDIR")
    workDir = os.path.realpath(workDirIni)
    for i in _workdirsMatix:
      os.environ["%s_WORKDIR" % i] = os.path.join(workDir, i)
    if workDirIni != workDir: #in case of links
      logger.warning("There is environment variable MATIX_WORKDIR as link:\n'%s'\n  ----> as real path:\n'%s'" % (workDirIni, workDir))
    else:
      logger.warning("Set environment variable MATIX_WORKDIR to:\n  '%s'" % workDir)
  else: 
    workDir = os.path.realpath(workDir)
  
  if not os.path.exists(workDir):
    os.makedirs(workDir)
    logger.warning('Create workdir ' + workDir)

  for i in _workdirsMatix:
    envName = "%s_WORKDIR" % i
    workDir = os.getenv(envName)
    if workDir is None:
      workDirByDefault = os.path.join(os.environ["DEFAULT_WORKDIR"] , i)
      logger.error('Environ variable unknown %s set default %s' % (envName, workDirByDefault))
      workDir = os.path.realpath(os.path.expandvars(workDirByDefault))
      os.environ[envName] = workDir


    workDirReal = os.path.realpath(os.path.expandvars(workDir))
    if not os.path.exists(workDirReal):
      os.makedirs(workDirReal)
      logger.warning('Create workdir\n%s\nas real path\n%s' % (workDir, workDirReal))

  _initializeWorkdirMatixDone = True
  return

def getHtmlDoc():
  """search for 'index.html' from relative salomepy directory path"""
  here = os.path.dirname(os.path.realpath(__file__))
  # PACKAGESPY salome
  path = os.path.realpath(os.path.join(here, "../../doc/index.html"))
  if os.path.isfile(path): return path
  # IRADINAGUI standalone
  path = os.path.realpath(os.path.join(here, "../doc/build/html/index.html"))
  if os.path.isfile(path): return path
  logger.error("Html doc as file .../index.html not found")
  return None

