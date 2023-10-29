#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
very elementary api for git repo commands
systematically calling subprocess.popen in current dir repo
user warning if stderr is not empty
user stuff to analyse stdout stderr "as it"

this is to avoid pip install gitpython or else in salome context (...for now)

example usage:
  import gitEmulate as git
  repo = git.Repo( '/tmp/repodir' )
  repo.isDir()
  repo.isGitRepoDir()
  repo._cmd("git init")
  stdout, stderr = repo._cmd("git status")
  etc...  
"""

import os
import sys
import subprocess as SP

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False

def root_path():
  return os.path.abspath(os.sep)

_forbidden = ["", None, root_path()] #risky for git init!!!!
_gitignoreDefault = """
*~
*.pyc
*.pyo
nohup.out

#created by highlight
highlight.css
"""

class Repo(object):
  def __init__(self, repo):
    """repo is a directory name, existing or not, with environ variable or not"""
    if repo in _forbidden: #avoid git init in "/"
      logger.critical("Problem path: '%s'" % repo)
      return
    self._envVarRepo = repo
    self._realRepo = self.interpretEnvVar(repo)
    if self._realRepo in _forbidden:
      raise Exception("Problem  path: '%s'" % self._realRepo)
    
  def _root_path():
    return os.path.abspath(os.sep)
    
  def _cmd(self, cmd, Verbose=verbose): 
    """
    all commands have Repo current directory
    returns (stdout, stderr)
    """
    if Verbose: print("\nRepo._cmd cmd = '%s'" % cmd)
    proc = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE, cwd=self._realRepo)
    stdout, stderr = proc.communicate()
    if Verbose: print("stdout: '%s'\nstderr: '%s'" % (stdout, stderr))
    return stdout, stderr

  def _cmdView(self, cmd, centralLogView, Verbose=False): 
    """
    all commands have Repo current directory
    display (stdout, stderr) in centralLogView by
    """
    if Verbose: print("\nRepo._cmdView cmd = '%s'" % cmd)
    centralLogView.launchCmdIntoPopen(cmd, cwd=self._realRepo)
    return "", "" #as _cmd returns #TODO get stdout, stderr

  def interpretEnvVar(self, aValueWithEnvVar):
    """
    resolve env variable as ${HOME} etc... 
    with os.path.expandvars
    """
    res = os.path.expandvars(aValueWithEnvVar)
    if "$" in res:
      logger.warning("Problem expandvars in: '%s' -> '%s'" % (aValueWithEnvVar, res))
    return res
    
  def getRealPath(self, aPathWithEnvVar):
    """
    resolve env variable as $HOME/toto etc... 
    with subprocess shell interpretation of env var
    """
    res = os.path.expandvars(aPathWithEnvVar)
    res = os.path.realpath(res)
    if verbose: print("getRealPath: '%s'->'%s'" % (aPathWithEnvVar, res))
    return res

  def isGitRepoDir(self):
    """
    if exist directory .git
    warning: may be sometimes creation multithread 'en cours'
    """
    if not os.path.isdir(self._realRepo): return False
    if not os.path.isdir(os.path.join(self._realRepo, ".git")): 
      #print "*****not a gitrepo******",os.path.join(self._realRepo, ".git")
      return False
    return True

  def isDir(self):
    if not os.path.isdir(self._realRepo): return False
    return True

  def createDir(self):
    if os.path.isdir(self._realRepo): 
      logger.warning("directory existing yet: '%s'" % self._realRepo)
      return True
    else:
      try: 
        os.makedirs(self._realRepo)
        return True
      except:
        logger.error("error creating directory: '%s'" % self._realRepo)
        return False

  def createGitignore(self, contents=None):
    name = os.path.join(self._realRepo, ".gitignore")
    try:
      if contents == None:
        conts = _gitignoreDefault
      else:
        conts = contents
      with open(name, "w") as f: f.write(conts)
    except:
      logger.error("error creating .gitignore file: '%s'" % name)
      return False
    return os.path.isfile(name)

    


