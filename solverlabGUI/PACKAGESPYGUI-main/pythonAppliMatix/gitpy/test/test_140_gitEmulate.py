#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import shutil
import unittest
import gitpy.gitEmulate as git

def root_path():
  return os.path.abspath(os.sep)

verbose = False

user = os.getenv("USER")

"""reading file from children directory ./tests_uranie"""
testDir = os.path.join(root_path() + "tmp" , user + "_gitEmulateTest")

#########################################
def toString(value):
  """
  assume conversion/decode if value is byte from unicode python2/3
  see https://docs.python.org/3/howto/unicode.html
  """
  typ = type(value)
  if typ is str:
    return value
  if typ is bytes:
    res = value.decode("utf-8", "ignore") # may be latin-1 better
    return res
  logger.warning("unexpected type %s for python2-3 string coding of ''" % (typ, value))
  return value

#########################################
class TestCase(unittest.TestCase):

  def test_005(self):
    if verbose: print("testDir %s" % testDir)
    if os.path.isdir(testDir): shutil.rmtree(testDir)
    if os.path.isfile(testDir): os.remove(testDir)
    self.assertFalse(os.path.exists(testDir))

  def test_010(self):
    if verbose: repo = git.Repo("") #get message
    if verbose: repo = git.Repo(None)
    if verbose: repo = git.Repo(root_path())

  def test_020(self):
    repo = git.Repo(testDir)
    self.assertFalse(repo.isDir())
    self.assertFalse(repo.isGitRepoDir())
    self.assertTrue(repo.createDir())
    if verbose: self.assertTrue(repo.createDir()) #get message
    self.assertTrue(repo.isDir())
    self.assertFalse(repo.isGitRepoDir())
    stdout, stderr = repo._cmd("git init")
    stdout = toString(stdout)
    # print("git init stdout", testDir, stdout)
    self.assertTrue(os.path.isdir(testDir+'/.git'))
    # stderr = toString(stderr) # get astuce: ... avoid stderr
    # self.assertEqual(stderr, "")
    # print("stdout:\n%s\nstderr:\n%s" % (stdout, stderr))
    # self.assertTrue("Initialized empty Git repository" in stdout or \
    #                 "Dépôt Git vide initialisé dans" in stdout)
    self.assertTrue(repo.createGitignore(contents = git._gitignoreDefault)) #default
    self.assertTrue(repo.isGitRepoDir())
    self.assertTrue(os.path.isdir(os.path.join(testDir, ".git")))

if __name__ == '__main__':
  verbose = False
  """
  True to see message
  CRITICAL : gitEmulate.py (69) : __init__ : Problem path: ''
  CRITICAL : gitEmulate.py (69) : __init__ : Problem path: 'None'
  CRITICAL : gitEmulate.py (69) : __init__ : Problem path: '/'
  WARNING : gitEmulate.py (136) : createDir : directory existing yet: '/tmp/wambeke_gitEmulateTest'
  """
  unittest.main()
  pass
