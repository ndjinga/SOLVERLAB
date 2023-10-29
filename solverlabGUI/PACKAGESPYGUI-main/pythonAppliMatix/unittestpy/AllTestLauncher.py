#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
Test all unittest files (test_*.py) existing in subdirectories of a root directory

| Example:
| >>  matix context
| >>  AllTestLauncher.py --rootPath='.' --pattern='test_???_*.py'
"""

"""
#### base of algorithm functionality

#see http://www.voidspace.org.uk/python/articles/introduction-to-unittest.shtml

import unittest

import test_something
import test_something2
import test_something3

loader = unittest.TestLoader()

suite = loader.loadTestsFromModule(test_something)
suite.addTests(loader.loadTestsFromModule(test_something2))
suite.addTests(loader.loadTestsFromModule(test_something3))

runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(suite)
"""


import os
import sys
import unittest
import traceback
from collections import OrderedDict
import argparse as AP

import glob
import fnmatch
import pprint as PP #pretty print

debug = False
verboseImport = True

# get path to origin sources
defaultdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

_user = os.getenv('USERNAME')
# wambeke is christian at home
_developers = ["christian", "wambeke",] #  ...

def errPrint(aStr):
  """stderr to avoid write in html or xml file log message"""
  sys.stderr.write(aStr + '\n')


if sys.path[0] != defaultdir:
  sys.path.insert(0, defaultdir)
  errPrint("WARNING  : sys.path prepend '%s'\n" % defaultdir)

try:
  # supposedly use xyz logger if defaultdir in PATH
  import xyzpy.loggingXyz as LOG
  logger = LOG.getLogger()
except:
  import logging as LOG
  logger = LOG.getLogger()
  logger.warning("use a default logger")

logger.warning("Python version is %s" % str(sys.version_info[0:3]))

try:
  import unittestpy.HTMLTestRunner as HTST
except:
  HTST = None
  logger.warning("""
no HTML output available.
try find 'test/unittestpy/HTMLTestRunner.py'
""")

try:
  import xmlrunner as XTST
except:
  XTST = None
  logger.warning("""\
no XML output available for unittest.
try 'pip install unittest-xml-reporting'
""")


###################################################################
def locate(pattern, root=os.curdir):
  """
  Locate all files matching supplied filename pattern in and below
  supplied root directory.
  """
  result = []
  logger.info("locate '%s' from '%s'\n" % (pattern, defaultdir))
  for path, dirs, files in os.walk(os.path.abspath(root)):
    for filename in fnmatch.filter(files, pattern):
      result.append( os.path.join(path, filename) )
  return result

def printEnv(search=""):
  """
  list all environment variables which name contains search string
  example: 
    import AllTestLauncher as ATL
    ATL.printEnv("ROOT_")
  """
  env=os.environ
  for i in sorted(env):
    if search in i:
      print(i) 

def grepInEnv(search=""):
  """
  list all environment variables which contains search string
  example: 
    import AllTestLauncher as ATL
    ATL.grepInEnv("XDATA")
  """
  env=os.environ
  for i in sorted(env):
     done=False
     for j in env[i].split(":"):
       if search in j:
           if not done:
             print(i+" contains ")
             done=True
           print("  "+j)
      
def format_exception(msg, limit=None, trace=None):
    """
    Format a stack trace and the exception information.
    as traceback.format_exception(),
    with all traceback only if user in ._developers
    """
    etype, value, tb = sys.exc_info()
    if _user in _developers:
      res = "\n" + msg
      if tb:
          res += "\nTraceback (most recent call last):\n"
          res += "".join(traceback.format_tb(tb, limit)) #[:-1])
      res += "\n<"
      res += "\n".join(traceback.format_exception_only(etype, value))
      return res
    else:
      res = "\n" + msg
      if tb:
          res += "\nTraceback:\n"
          res += "".join(traceback.format_tb(tb, limit)[-1:]) #[:-1])
      res += "\n<"
      res += "".join(traceback.format_exception_only(etype, value))
      return res

###################################################################
def runOnArgs(args):
  """
  launch tests on args.pattern files
  """
  import logging as LOGI
  logger = LOGI.getLogger()
  logger.trace = logger.debug
  logger.step = logger.info

  fromFileOrPath = args.rootPath
  fileTestPattern = args.pattern
  if fromFileOrPath == None:
    directory, name = os.path.split( os.path.realpath( __file__ ) )
  else:
    if os.path.isdir(fromFileOrPath):
      directory, name = (fromFileOrPath, None)
      fileTestPatternCurrent = fileTestPattern
    elif os.path.isfile(fromFileOrPath):
      directory, name = os.path.split( os.path.realpath( fromFileOrPath ) )
      fileTestPatternCurrent = name
    else:
      mess = "Cannot get file or directory '%s'" % fromFileOrPath
      errPrint("ERROR : " + mess)
      return None
      #raise Exception("Cannot get file or directory '%s'" % fromFileOrPath)

  #files = glob.glob(os.path.join(directory, "*Test.py"))
  files = locate(fileTestPatternCurrent, directory)
  files = sorted(files, key=lambda bas: os.path.basename(bas))

  filesForTest = OrderedDict()

  if verboseImport:
    errPrint("\n***************** try import tagged file(s) '%s' ***************** \n" % \
             fileTestPattern)
  for aFile in files:
    aDir, aName = os.path.split(aFile)
    aImport, ext = os.path.splitext(aName)
    
    try:
      if aFile in list(filesForTest.keys()):
        print("WARNING: imported yet: %s" % aFile)
      else:
        sys.path.insert(0, aDir) # beware of sys.path appends at each file (for robustness)
        done = True
        if verboseImport: errPrint("try import '%s'" % aImport)
        aModule = __import__(aImport, globals(), locals(), []) 
        del sys.path[0]
        done = False
        filesForTest[aFile] = (aImport, aModule)
    except Exception as e:
      if done: 
        del sys.path[0] # beware of sys.path appends at each file (for robustness)
        done = False
      msg = "ERROR : AllTestLauncher: import %s:" % aFile
      err = format_exception(msg)
      errPrint(err)
      continue

  listfilesForTest = list(filesForTest.keys()) # avoid type odict_keys python3
  result = None

  errPrint("\n\n***************** AllTestLauncher test files(s) *****************\n\n%s\n" % \
           PP.pformat(listfilesForTest) )
  
  if len(listfilesForTest) == 0: 
    if debug: errPrint("WARNING : AllTestLauncher: empty list of test files")
    return None

  loader = unittest.TestLoader()
  suite = None

  for i,k in enumerate(listfilesForTest):
    if debug: errPrint("Test : %s %s" % (i, k))
    if i == 0:
      suite = loader.loadTestsFromModule( filesForTest[k][1] )
      pass
    else:
      suite.addTests( loader.loadTestsFromModule( filesForTest[k][1] ) )
      pass

  if args.type == "std": 
    runner = unittest.TextTestRunner(verbosity=args.verbosity)
  elif args.type == "html": 
    runner = HTST.HTMLTestRunner(verbosity=args.verbosity, )
  elif args.type == "xml": 
    if args.name == 'stdout':
      #all-in-one xml output at 'sys.stdout' for pipe redirection
      runner = XTST.XMLTestRunner(verbosity=args.verbosity, output=sys.stdout)
    else:
      #one file xml per test in suite in args.name directory
      runner = XTST.XMLTestRunner(verbosity=args.verbosity, output=args.name)
  else:
    errPrint("ERROR : unknown type of output: '%s'" % args.type)
    return None    
    
  if suite != None: result = runner.run(suite)
  return result

###################################################################
def runFromEnvVar(envVar, fileTestPattern="test_*.py"):
  """
  example: 
    import AllTestLauncher as ATL
    ATL.runFromEnvVar("MICROGEN_ROOT_DIR")
    ATL.runFromEnvVar("MICROGEN_ROOT_DIR", "aggregate_*GJKTest.py")
  """
  env=os.environ
  res = []
  for i in sorted(env):
    if envVar in i:
      res.append(i)
  if len(res) > 1:
    mess = "multiple environment variable for '%s': %s" % (envVar, str(res))
    errPrint("ERROR : " + mess)
    return None
  if len(res) < 1:
    mess = "no environment variable for '%s'" % (envVar)
    errPrint("ERROR : " + mess)
    return None
  res = res[0]
  tmp = env[res].split(":")
  if len(tmp) > 1:
    mess = "need only one path in environment variable '%s'" % (res)
    errPrint("ERROR : " + mess)
    return None  
  run(fromFileOrPath=env[res], fileTestPattern=fileTestPattern)


###################################################################
def getParser():
  parser = AP.ArgumentParser(description='launch All python tests', argument_default=None)

  parser.add_argument(
    '-d', '--debug', 
    help='set debug mode, more verbose',
    action='store_true',
  )
  parser.add_argument(
    '-v', '--verbosity', 
    help='set verbosity of unittests [0|1|2...]',
    default=2,
    metavar='int'
  )
  parser.add_argument(
    '-r', '--rootPath', 
    help="""\
dir name with absolute or relative path stand for root directory
of recursive searching unittest python files
""",
   default=defaultdir,
   metavar='dirPath'
  )
  parser.add_argument(
    '-p', '--pattern', 
    help="file pattern for unittest files ['test_*.py'|'*Test.py'...]",
    default="test_???_*.py", # as alphabetical ordered test site
    metavar='filePattern'
  )
  parser.add_argument(
    '-t', '--type', 
    help="type of output: ['std'(standart ascii)|'xml'|'html']",
    default="std",
    choices=['std', 'xml', 'html'],
    metavar='outputType'
  )
  parser.add_argument(
    '-n', '--name', 
    help="""\
(only for type xml)
name of directory output: ['test_reports'|...].
If name = 'stdout' then all-in-one xml output at 'sys.stdout'. For pipe redirection:
'>> AllTestLauncher.py -t xml -n stdout > tmp.xml'
""",
    default="test_reports",
    metavar='dirName'
  )
  return parser

#export PATH=defaultdir:${PATH}

###################################################################
if __name__ == '__main__':
  # Make the src & command package accessible from all code
  # as export PYTHONPATH=defaultdir:${PYTHONPATH}
  # https://docs.python.org/2/library/os.html
  # On some platforms, including FreeBSD and Mac OS X, 
  # setting environ may cause memory leak
  # so use sys.path

  # errPrint("INFO : AllTestLauncher sys.path:\n'%s'" % PP.pformat(sys.path)

  args = getParser().parse_args(sys.argv[1:])
  debug = args.debug
  if debug:
    print("INFO: args:\n  %s" % PP.pformat(args))
  directory = os.path.realpath(args.rootPath) # supposed to be root of a package
  if directory not in sys.path: # to assume import from directory are ok
    sys.path.insert(0, directory)
    print("INFO : AllTestLauncher sys.path.insert %s" % directory)

  if debug:
    print("DEBUG : AllTestLauncher environ:\n%s" % PP.pformat(dict(os.environ)))
  runOnArgs(args)
  # useless del sys.path[0] if done



