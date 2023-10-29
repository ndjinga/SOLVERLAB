#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
Code during executing thread
From the Open Group page for the open() function
O_NONBLOCK
  When opening a FIFO with O_RDONLY or O_WRONLY set:
  If O_NONBLOCK is set:
      An open() for reading only will return without delay. An open()
      for writing only will return an error if no process currently
      has the file open for reading.

  If O_NONBLOCK is clear:
      An open() for reading only will block the calling thread until a
      thread opens the file for writing. An open() for writing only
      will block the calling thread until a thread opens the file for
      reading.
"""

import os
import sys
import stat
import errno
import time
import pprint as PP

import subprocess as SP
import threading as THRD

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False

_aPipeNameDefault = "pipefifo_001" # default choice for a named pipe file


#############################################
def os_write(aPipe, aString):
  """
  assume os.write without error
  'a bytes-like object is required, not 'str'' in python 3
  """
  _python = sys.version_info[0] # python 2 or 3
  # print("python %s" % _python)
  if _python == 2:
    os.write(aPipe, aString)
    return
  else: # python 3 needs bytes
    aBytes = bytes(aString, "utf_8")
    os.write(aPipe, aBytes)
    return

#############################################
def toString(strOrBytes):
  """
  assume os.read as 'str', not 'bytes-like object' in python 3
  """
  if type(strOrBytes) is bytes:
    return strOrBytes.decode("utf-8")
  else:
    return strOrBytes


_testWriteBlocking = r'''#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this file is for test
automatically created by
@FILE@

to launch test:
  @LAUNCHCMD@
"""

import os
import sys
import time

verbose = True

#############################################
def os_write(aPipe, aString):
  """
  assume os.write without error
  'a bytes-like object is required, not 'str'' in python 3
  """
  _python = sys.version_info[0] # python 2 or 3
  # print("python %s" % _python)
  if _python == 2:
    os.write(aPipe, aString)
    return
  else: # python 3 needs bytes
    aBytes = bytes(aString, "utf_8")
    os.write(aPipe, aBytes)
    return

#############################################
# write pipe
aFile = "@PIPE@"
aTime = time.strftime("%Hh %Mm %Ss")
mess = "Hello from _testPipeFifo at %s\n" % aTime
# os open
pipeout = os.open(aFile, os.O_WRONLY)   # blocking if not open "r" at other pipe side
if verbose: print("_testPipeFifo open pipe write done")
os_write(pipeout, mess)
if verbose: print("_testPipeFifo write pipe send: '%s'" % mess[:-1])
os.close(pipeout)
print('_testPipeFifo end write')

'''

_testReadBlocking = r'''#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this file is for test
automatically created by
@FILE@
"""

import os
import sys

verbose = True

#############################################
def toString(strOrBytes):
  """
  assume os.read as 'str', not 'bytes-like object' in python 3
  """
  if type(strOrBytes) is bytes:
    return strOrBytes.decode("utf-8")
  else:
    return strOrBytes


#############################################
# read pipe
aFile = "@PIPE@"
# python open
pipein = open(aFile, 'r')  # blocking if not open "w" at other pipe side
if verbose: print('_testPipeFifo open pipe read done')
line = pipein.readline()   # blocks until data sent
if verbose: print("_testPipeFifo read pipe receive: '%s'" % line[:-1])
pipein.close()
print('_testPipeFifo end read')

'''

_testWriteNoBlocking = r'''#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this file is for test
automatically created by
@FILE@

to launch test:
  @LAUNCHCMD@
"""

import os
import sys
import time
import errno

warning = True

#############################################
def os_write(aPipe, aString):
  """
  assume os.write without error
  'a bytes-like object is required, not 'str'' in python 3
  """
  _python = sys.version_info[0] # python 2 or 3
  # print("python %s" % _python)
  if _python == 2:
    os.write(aPipe, aString)
    return
  else: # python 3 needs bytes
    aBytes = bytes(aString, "utf_8")
    os.write(aPipe, aBytes)
    return

#############################################
# write pipe
aFile = "@PIPE@"
# os open
try:
  pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)   # no blocking if not open "r" at other pipe side
except OSError as exc: # seems to be when pipe is not open read
  ermsg = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno]) + \
          "\n  _testPipeFifo open pipe write needs to be somewhere previously open read, (fixed for you):\n  %s" % aFile
  print(ermsg)
  pipein = os.open(aFile, os.O_RDONLY | os.O_NONBLOCK)   # open named pipe fifo
  pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)  # no problem if not open "r" at other pipe side

if warning: print("_testPipeFifo open pipe write done")

# first send
aTime = time.strftime("%Hh %Mm %Ss")
mess = "Hello 1 from _testPipeFifo at %s\n" % aTime
os_write(pipeout, mess)
if warning: print("_testPipeFifo write pipe send: '%s'" % mess[:-1])

# second send
time.sleep(4)
aTime = time.strftime("%Hh %Mm %Ss")
mess = "Hello 2 from _testPipeFifo at %s\n" % aTime
os_write(pipeout, mess)
if warning: print("_testPipeFifo write pipe send: '%s'" % mess[:-1])

os.close(pipeout)

# third send
time.sleep(3)
# os open again
pipeout = os.open(aFile, os.O_WRONLY)   # blocking if not open "r" at other pipe side
if warning: print("_testPipeFifo open again pipe write done")
aTime = time.strftime("%Hh %Mm %Ss")
mess = "Hello 3 from _testPipeFifo at %s\n" % aTime
os_write(pipeout, mess)
if warning: print("_testPipeFifo write pipe send: '%s'" % mess[:-1])
os.close(pipeout)

print('_testPipeFifo end write NONBLOCK')

'''

_testWriteNoBlockingSimple = r'''#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this file is for test
automatically created by
@FILE@

to launch test:
  @LAUNCHCMD@
"""

import os
import sys
import time
import errno

warning = True

#############################################
def os_write(aPipe, aString):
  """
  assume os.write without error
  'a bytes-like object is required, not 'str'' in python 3
  """
  _python = sys.version_info[0] # python 2 or 3
  # print("python %s" % _python)
  if _python == 2:
    os.write(aPipe, aString)
    return
  else: # python 3 needs bytes
    aBytes = bytes(aString, "utf_8")
    os.write(aPipe, aBytes)
    return

#############################################
# write pipe
aFile = "@PIPE@"
# os open
try:
  pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)   # no blocking if not open "r" at other pipe side
except OSError as exc: # seems to be when pipe is not open read
  ermsg = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno]) + \
          "\n  _testPipeFifo open pipe write needs to be somewhere previously open read, (fixed for you):\n  %s" % aFile
  print(ermsg)
  pipein = os.open(aFile, os.O_RDONLY | os.O_NONBLOCK)   # open named pipe fifo
  pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)  # no problem if not open "r" at other pipe side

if warning: print("_testPipeFifo open pipe write done")

# first send
aTime = time.strftime("%Hh %Mm %Ss")
mess = "# Hello 1 from _testPipeFifo at %s\n" % aTime
os_write(pipeout, mess)
if warning: print("_testPipeFifo write pipe send: '%s'" % mess[:-1])
os.close(pipeout)

print('_testPipeFifo end write NONBLOCK simple')

'''



_testReadNoBlocking = r'''#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this file is for test
automatically created by
@FILE@
"""

import os
import sys
import errno
import time

verbose = True

#############################################
def toString(strOrBytes):
  """
  assume os.read as 'str', not 'bytes-like object' in python 3
  """
  if type(strOrBytes) is bytes:
    return strOrBytes.decode("utf-8")
  else:
    return strOrBytes

#############################################
# read pipe
aFile = "@PIPE@"

pipein = os.open(aFile, os.O_RDONLY | os.O_NONBLOCK) # open fifo pipe file as fd
if verbose: print('_testPipeFifo open pipe read NONBLOCK done')

i = 0
imax = 15 #as seconds total loop
while True:
  i += 1
  if i > imax: #max loop
   print("_testPipeFifo read NONBLOCK is stopping as expected at %i" % imax)
   break

  try:
    #as 1024 bytes per second max with time.sleep(1)
    line = os.read(pipein, 1024)
    line = toString(line)
  except OSError as exc: # seems to be when write is open but empty
    ermsg = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno])
    print(ermsg)
    line = ""
    pass

  if line == "": # seems to be when write is closed
    print("_testPipeFifo read pipe receive: 'NoMessage' i=%03d " %i)
  else:
    print("_testPipeFifo read pipe receive: '%s' i=%03d" % (line, i))
  time.sleep(1)

os.close(pipein)
print('_testPipeFifo end read NONBLOCK')

'''

_testWriteRootNoBlocking = r'''
{
  /*
  this file is for test
  automatically created by
  @FILE@

  to launch test:
  @LAUNCHCMD@
  */

  TDatime theTime;
  ofstream sortie;

  sortie.open("@PIPE@", ios::out);

  theTime.Set();
  sortie << "Hello 1 from _testPipeFifo at " << theTime.AsString() << endl;

  gSystem->Sleep(2000); // wait 2 sec

  theTime.Set();
  sortie << "Hello 2 from _testPipeFifo at " << theTime.AsString() << endl;

  sortie.close();
  cout << "_testPipeFifo end write NONBLOCK" << endl;
}

'''


#http://www.cplusplus.com/forum/unices/107729/
_testWriteMessageRootNoBlocking = r'''
{
  /*
  this file is for test
  automatically created by
  @FILE@

  to launch test:
  @LAUNCHCMD@
  */

  TDatime theTime;
  ofstream sortie;

  sortie.open("@PIPE@", ios::out);

  theTime.Set();
  sortie << "@MESSAGE@ # from _testPipeFifo at " << theTime.AsString() << endl;
  gSystem->Sleep(1000); // wait 2 sec

  theTime.Set();
  sortie << "@MESSAGE@ # from _testPipeFifo at " << theTime.AsString() << endl;
  gSystem->Sleep(1000); // wait 2 sec

  sortie.close();
  cout << "_testPipeFifo end write NONBLOCK" << endl;
}

'''

_old____testWriteMessageRootNoBlocking = r'''
{
  /*
  this file is for test
  automatically created by
  @FILE@

  to launch test:
  @LAUNCHCMD@
  */

  TString mess = "@MESSAGE@";
  TDatime theTime;
  Int_t ok = 0;

  printf("  Trying to open.\n");
  FILE *fd = fopen("@PIPE@", "w");
  printf("  Connection succeeded.\n");

  theTime.Set();
  // mandatory end line in mess for readline()
  mess = mess + " # from _testPipeFifo at " + theTime.AsString() + "\n";
  ok = fwrite( mess.Data(), 1, mess.Length(), fd );
  cout << "  write " << mess.Length() << " bytes, cr=" << ok << ":" << endl;
  cout << mess << endl;

  // gSystem->Sleep(2000); // wait 2 sec

  cerr << "before close" << endl;
  ok = fclose(fd);
  // ok =-1 as stderr *** Break *** write on a pipe with no one to read it
  printf("  _testPipeFifo end write NONBLOCK FILE.\n");
  cerr << "end prog ok=" << ok << endl;
}

'''

_Problem1_____testWriteMessageRootNoBlocking = r'''
{
  /*
  this file is for test
  automatically created by
  @FILE@

  to launch test:
  @LAUNCHCMD@
  */

  TDatime theTime;
  ofstream sortie;

  sortie.open("@PIPE@", ios::out);

  theTime.Set();
  sortie << "@MESSAGE@ # from _testPipeFifo at " << theTime.AsString() << endl;

  gSystem->Sleep(2000); // wait 2 sec

  sortie.close();
  cout << "_testPipeFifo end write NONBLOCK" << endl;
}

'''

_Problem2_____testWriteMessageRootNoBlocking = r'''
{
  /*
  this file is for test
  automatically created by
  @FILE@

  to launch test:
  @LAUNCHCMD@
  */

  TDatime theTime;
  int pipe;

  printf("Trying to open.\n");
	if ( pipe = open("@PIPE@", O_WRONLY ) < 1 ) {
		printf("Error: open()");
		return 1;
	}

  printf("Connection succeeded.\n");

  theTime.Set();
  sortie << "@MESSAGE@ # from _testPipeFifo at " << theTime.AsString() << endl;
  if ( write( pipe, &number, sizeof( int ) ) < 0 ) {
		printf(("Error: write()");
		return 1;
	 }

  gSystem->Sleep(2000); // wait 2 sec

  sortie.close();
  cout << "_testPipeFifo end write NONBLOCK" << endl;
}
'''

_testReadRootNoBlocking = '''

TODO _testReadRootNoBlocking

'''

###########################################################################
def isFilePipeFifo(realpath):
  """return False if inexisting or existing as classical file"""
  try:
    if os.path.exists(realpath): # check if file is a named pipe (fifo)
      return stat.S_ISFIFO(os.stat(realpath).st_mode)
  except: #case realpath=None or else...
    pass
  return False

def chmodax(aFile):
  """chmod a+x path"""
  st = os.stat(aFile)
  os.chmod(aFile, st.st_mode | 0o111)


###########################################################################
class PipeFifo(object):
  def __init__(self, nameFifo=None):
    self._nameFifo = None #private
    if nameFifo != None:
      self.setNameFifo(nameFifo)
    self.pipein = None #to get writeNoBlocking pipe needs to be previously open read

  def getNameFifo(self):
    return self._nameFifo

  def getBaseNameFifo(self):
    try:
      _, base = os.path.split(self.getNameFifo())
    except:
      base = None
    return base

  def setNameFifo(self, nameFifo):
    """only one set"""
    if self._nameFifo != None:
      raise Exception("PipeFifo: pipe set yet:\n%s" % self)
    if not os.path.isfile(nameFifo):
      self._createPipeFifo(nameFifo)

  def _createPipeFifo(self, nameFifo, createDir=True):
    """only one set for one instance"""
    realpath = os.path.realpath(nameFifo)
    aDir, aFile = os.path.split(realpath)
    if not createDir:
      if not os.path.isdir(aDir):
        raise Exception("PipeFifo: directory not existing: %s" % aDir)
    else:
      if not os.path.isdir(aDir):
        os.makedirs(aDir)
        if verbose: logger.warning("create directory for named pipe fifo: %s" % aDir)
    if os.path.exists(realpath): # check if file is a named pipe (fifo)
      if not stat.S_ISFIFO(os.stat(realpath).st_mode):
        raise Exception("PipeFifo: existing file is not named pipe fifo: %s" % realpath)
    else:
      os.mkfifo(realpath)        # create the named pipe file
      if verbose: logger.info("create named pipe fifo %s" % realpath)
    self._nameFifo = realpath
    time.sleep(.1) #cool attitude precaution

  def isPipeFifo(self):
    """return False in case of external process of file copy or recreate"""
    return isFilePipeFifo(self.getNameFifo())

  def __repr__(self):
    return "%s(nameFifo='%s', isPipeFifo='%s')" % \
           ( self.__class__.__name__, self.getNameFifo(), self.isPipeFifo() )

  ########################################################################################
  # WARNING:
  # user have to use these read or write methods exclusively in
  # distinct threads or distinct subprocess
  ########################################################################################
  def writeBlocking(self, mess):
    """could be blocking"""
    aFile = self.getNameFifo()
    if self.isPipeFifo():
      # os open
      if verbose: logger.info('PipeFifo.writeBlocking(): try open pipe write')
      pipeout = os.open(aFile, os.O_WRONLY)   # blocking if not open "r" at other pipe side
      if verbose: logger.info("PipeFifo.writeBlocking(): open pipe write done")
      os_write(pipeout, mess)
      if verbose: logger.info("PipeFifo.writeBlocking(): write pipe send: '%s'" % mess[:-1])
      os.close(pipeout)
      return True
    else:
      raise Exception("PipeFifo.writeBlocking(): file is not named pipe fifo: %s" % aFile)

  def writeNoBlocking(self, mess, Verbose=True):
    """could be blocking"""
    aFile = self.getNameFifo()
    if self.isPipeFifo():
      # os open
      if verbose: logger.info('PipeFifo.writeNoBlocking(): try open pipe write')
      try:
        pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)   # no blocking if not open "r" at other pipe side
      except OSError as exc: # seems to be when pipe is not open read
        if Verbose:
          ermsg = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno]) + \
                  "\n  pipe needs to be somewhere previously open read, (fixed for you):\n  %s" % aFile
          logger.warning(ermsg)
        self.pipein = os.open(aFile, os.O_RDONLY | os.O_NONBLOCK)   # open named pipe fifo
        pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)  # no problem if not open "r" at other pipe side
      if verbose: logger.info("PipeFifo.writeNoBlocking(): open pipe write done")
      os_write(pipeout, mess)
      if verbose: logger.info("PipeFifo.writeNoBlocking(): write pipe send: '%s'" % mess)
      os.close(pipeout)
      #if self.pipein != None: os.close(self.pipein) #doing that seems clear pipe
      return True
    else:
      raise Exception("PipeFifo.writeNoBlocking(): file is not named pipe fifo: %s" % aFile)

  def readlineBlocking(self):
    """could be blocking"""
    aFile = self.getNameFifo()
    if self.isPipeFifo():
      # python open
      if verbose: logger.info('PipeFifo.readlineBlocking(): try open pipe read %s' % aFile)
      pipein = open(aFile, 'r')  # blocking if not open "w" at other pipe side
      if verbose: logger.info('PipeFifo.readlineBlocking(): open pipe read done')
      line = pipein.readline()   # blocks until data sent
      if verbose: logger.info("PipeFifo.readlineBlocking(): pipe receive: '%s'" % line[:-1])
      pipein.close()
      return toString(line)
    else:
      raise Exception("PipeFifo.readlineBlocking(): file is not named pipe fifo: %s" % aFile)

  def readNoBlocking(self, buflen=1024):
    """read default max buflen=1024 bytes at each call"""
    aFile = self.getNameFifo()
    if self.isPipeFifo():
      if verbose: logger.info('PipeFifo.readNoBlocking(): try open pipe read')
      pipein = os.open(aFile, os.O_RDONLY | os.O_NONBLOCK) # open named pipe fifo
      if verbose: logger.info('PipeFifo.readNoBlocking(): open pipe read done')
      try:
        #as buflen bytes max read
        line = os.read(pipein, buflen)
      except OSError as exc: # seems to be when write is open but empty
        #line = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno])
        line = ""
        pass

      #seems to be line="" when write is closed
      if verbose: logger.info("PipeFifo.readNoBlocking(): pipe receive: '%s'" % line)
      os.close(pipein)
      return toString(line)
    else:
      raise Exception("PipeFifo.readNoBlocking(): file is not named pipe fifo: %s" % aFile)

  def readBlocking(self, buflen=1024):
    """read default max buflen=1024 bytes at each call"""
    aFile = self.getNameFifo()
    if self.isPipeFifo():
      if verbose: logger.info('PipeFifo.readBlocking(): try open pipe read')
      pipein = os.open(aFile, os.O_RDONLY) # open named pipe fifo
      if verbose: logger.info('PipeFifo.readBlocking(): open pipe read done')
      try:
        #as buflen bytes max read
        line = os.read(pipein, buflen)
      except OSError as exc: # seems to be when write is open but empty
        #line = "OSError %s (errno.%s)" % (exc.errno, errno.errorcode[exc.errno])
        line = ""
        pass

      #seems to be line="" when write is closed
      if verbose: logger.info("PipeFifo.readBlocking(): pipe receive: '%s'" % line)
      os.close(pipein)
      return toString(line)
    else:
      raise Exception("PipeFifo.readNoBlocking(): file is not named pipe fifo: %s" % aFile)


  ##################################################################################
  # tests methods (only for example and use unittest)
  ##################################################################################
  def _testPipeFifoBlocking(self):
    """short test of fifo, blocking mode, no thread, by subprocess"""
    aDir, aFile = os.path.split(self.getNameFifo())
    fileWrite = os.path.join(aDir, "testPipeFifoWrite.py")
    fileRead = os.path.join(aDir, "testPipeFifoRead.py")

    cmd = r"""
  # background first launch because interprocess blocking open/read/write test
  cd %s
  ls -alt
  %s &
  %s
""" % (aDir, fileWrite, fileRead)

    testWrite = _testWriteBlocking
    testWrite = testWrite.replace("@FILE@", os.path.realpath(__file__))
    testWrite = testWrite.replace("@PIPE@", self.getNameFifo())
    testWrite = testWrite.replace("@LAUNCHCMD@", cmd)

    testRead = _testReadBlocking
    testRead = testRead.replace("@FILE@", os.path.realpath(__file__))
    testRead = testRead.replace("@PIPE@", self.getNameFifo())

    with open(fileWrite, 'w') as f: f.write(testWrite)
    with open(fileRead, 'w') as f: f.write(testRead)
    chmodax(fileWrite)
    chmodax(fileRead)

    if verbose:
      logger.info('file Write: %s\n%s' % (fileWrite, testWrite))
      logger.infot('file Read: %s\n%s' % (fileRead, testRead))
      logger.info("launch ipc blocking simple test:%s" % cmd)

    stdout, stderr = SP.Popen(cmd, shell=True, cwd=aDir, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
    stdout = toString(stdout)
    stderr = toString(stderr)
    if stderr == "": stderr = "None"
    if verbose: logger.info("_testPipeFifoBlocking result:\nstdout: %s\nstderr: %s" % (stdout, stderr))
    if "read pipe receive: 'Hello" in stdout and \
       "write pipe send: 'Hello" in stdout:
      return True
    else:
      return False

  def _testPipeFifoNoBlocking(self):
    """short test of fifo, No blocking mode, no thread, by subprocess"""
    aDir, aFile = os.path.split(self.getNameFifo())
    fileWrite = os.path.join(aDir, "testPipeFifoWriteNoBlocking.py")
    fileRead = os.path.join(aDir, "testPipeFifoReadNoBlocking.py")

    cmd = r"""
  # background first launch because interprocess blocking open/read/write test
  cd %s
  ls -alt
  %s &
  %s
""" % (aDir, fileWrite, fileRead)

    testWrite = _testWriteNoBlocking
    testWrite = testWrite.replace("@FILE@", os.path.realpath(__file__))
    testWrite = testWrite.replace("@PIPE@", self.getNameFifo())
    testWrite = testWrite.replace("@LAUNCHCMD@", cmd)

    testRead = _testReadNoBlocking
    testRead = testRead.replace("@FILE@", os.path.realpath(__file__))
    testRead = testRead.replace("@PIPE@", self.getNameFifo())

    with open(fileWrite, 'w') as f: f.write(testWrite)
    with open(fileRead, 'w') as f: f.write(testRead)
    chmodax(fileWrite)
    chmodax(fileRead)

    if verbose:
      logger.info('file Write: %s\n%s' % (fileWrite, testWrite))
      logger.info('file Read: %s\n%s' % (fileRead, testRead))
      logger.info("launch ipc NO blocking simple test:%s" % cmd)

    stdout, stderr = SP.Popen(cmd, shell=True, cwd=aDir, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
    stdout = toString(stdout)
    stderr = toString(stderr)
    if stderr == "": stderr = "None"
    if verbose:
      logger.info("_testPipeFifoNoBlocking result:\nstdout: %s\nstderr: %s" % (stdout, stderr))
    if "read pipe receive: 'Hello 1" in stdout and \
       "read pipe receive: 'Hello 2" in stdout and \
       "end read NONBLOCK" in stdout and \
       "end write NONBLOCK" in stdout:
      return True
    else:
      return False

  def _testPipeFifoOnlyWriteNoBlocking(self, message='# None', Verbose=True):
    """short test of fifo only write, No blocking mode, no thread, by subprocess"""
    aDir, aFile = os.path.split(self.getNameFifo())
    fileWrite = os.path.join(aDir, "testPipeFifoWriteNoBlocking.py")

    cmd = r"""
  # background first launch because interprocess blocking open/read/write test
  cd %s
  ls -alt
  %s
""" % (aDir, fileWrite)

    testWrite = _testWriteNoBlockingSimple
    testWrite = testWrite.replace("@FILE@", os.path.realpath(__file__))
    testWrite = testWrite.replace("@PIPE@", self.getNameFifo())
    testWrite = testWrite.replace("@LAUNCHCMD@", cmd)

    with open(fileWrite, 'w') as f: f.write(testWrite)
    chmodax(fileWrite)

    if Verbose:
      #print('file Write: %s\n%s' % (fileWrite, testWrite))
      logger.info("launch ipc NO blocking simple test:%s" % cmd)

    stdout, stderr = SP.Popen(cmd, shell=True, cwd=aDir, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
    stdout = toString(stdout)
    stderr = toString(stderr)
    if stderr == "": stderr = "None"
    if Verbose:
      logger.info("_testPipeFifoOnlyWriteNoBlocking result:\nstdout: %s\nstderr: %s" % (stdout, stderr))
    if "end write NONBLOCK" in stdout:
      return True
    else:
      return False

  def _testPipeFifoRootNoBlocking(self, Verbose=False):
    """short test of fifo, No blocking mode, no thread, by subprocess"""
    aDir, aFile = os.path.split(self.getNameFifo())
    fileWrite = os.path.join(aDir, "testPipeFifoWriteRootNoBlocking.C")
    fileRead = os.path.join(aDir, "testPipeFifoReadNoBlocking.py")

    cmd = r"""
  # background first launch because interprocess blocking open/read/write test
  cd %s
  ls -alt
  which root
  root -l -q %s &
  %s
""" % (aDir, fileWrite, fileRead)

    testWrite = _testWriteRootNoBlocking
    testWrite = testWrite.replace("@FILE@", os.path.realpath(__file__))
    testWrite = testWrite.replace("@PIPE@", self.getNameFifo())
    testWrite = testWrite.replace("@LAUNCHCMD@", cmd)

    testRead = _testReadNoBlocking # _testReadRootNoBlocking
    testRead = testRead.replace("@FILE@", os.path.realpath(__file__))
    testRead = testRead.replace("@PIPE@", self.getNameFifo())

    with open(fileWrite, 'w') as f: f.write(testWrite)
    with open(fileRead, 'w') as f: f.write(testRead)
    chmodax(fileWrite)
    chmodax(fileRead)

    if verbose:
      logger.info('file Write: %s\n%s' % (fileWrite, testWrite))
      logger.info('file Read: %s\n%s' % (fileRead, testRead))
      logger.info("launch ipc Root NO blocking simple test:%s" % cmd)

    stdout, stderr = SP.Popen(cmd, shell=True, cwd=aDir, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
    stdout = toString(stdout)
    stderr = toString(stderr)
    if stderr == "": stderr = "None"
    if Verbose:
      logger.info("_testPipeFifoRootNoBlocking result:\nstdout: %s\nstderr: %s" % (stdout, stderr))
    if "read pipe receive: 'Hello 1" in stdout and \
       "read pipe receive: 'Hello 2" in stdout and \
       "end read NONBLOCK" in stdout and \
       "end write NONBLOCK" in stdout:
      return True
    else:
      return False

  def _testPipeFifoOnlyWriteRootNoBlocking(self, message='# None', Verbose=False):
    """
    short test of fifo, No blocking mode, no thread, by subprocess
    default message '# None' python compliant
    """
    aDir, aFile = os.path.split(self.getNameFifo())
    fileWrite = os.path.join(aDir, "testPipeFifoWriteRootNoBlocking.C")

    cmd = r"""
  # first launch interprocess no blocking write test
  cd %s
  pwd
  # ls -alt
  # which root
  afile=%s
  # cat $afile
  root -l -q $afile & # root write pipe
  echo 'WARNING: this script needs somebody reading pipe'
""" % (aDir, fileWrite)

    testWrite = _testWriteMessageRootNoBlocking
    testWrite = testWrite.replace("@FILE@", os.path.realpath(__file__))
    testWrite = testWrite.replace("@PIPE@", self.getNameFifo())
    testWrite = testWrite.replace("@LAUNCHCMD@", cmd)
    testWrite = testWrite.replace("@MESSAGE@", message)

    with open(fileWrite, 'w') as f: f.write(testWrite)
    chmodax(fileWrite)

    if Verbose:
      #print('file Write: %s\n%s' % (fileWrite, testWrite))
      logger.info('root write message: "%s"' % message)
      logger.info("launch ipc Root write NO blocking simple test:%s" % cmd)

    stdout, stderr = SP.Popen(cmd, shell=True, cwd=aDir, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
    stdout = toString(stdout)
    stderr = toString(stderr)

    if stderr == "": stderr = "None"
    if Verbose:
      logger.info("_testPipeFifoOnlyWriteRootNoBlocking result:\nstdout: %s\nstderr: %s" % (stdout, stderr))

    if "WARNING: this script needs somebody reading pipe" in stdout and \
       "end write NONBLOCK" in stdout:
      return True
    else:
      return False

#######################################################
#
# Threading named pipe fifo
#
# Discussions criticizing Python often talk about how it is difficult
# to use Python for multithreaded work, pointing fingers at what
# is known as the global interpreter lock (affectionately referred to as the “GIL”)
# that prevents multiple threads of Python code from running simultaneously.
# Due to this, the threading module doesn’t quite behave the way
# you would expect it to if you’re not a Python developer and
# you are coming from other languages such as C++ or Java
# see:
# https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python
#
#######################################################

def getListOfThreadsAlive():
  return THRD.enumerate()

def getStrListOfThreadsAlive():
  return "listOfThreadsAlive:\n%s" % PP.pformat(THRD.enumerate())

def printListOfThreadsAlive():
  logger.info(getStrListOfThreadsAlive())

#######################################################
class ThreadPipeFifoBase(THRD.Thread):
  """
  base class for
  Thread emit signal when send/receive as write/read string from named pipe fifo
  """
  noThread = [1]

  def __init__(self):
    super(ThreadPipeFifoBase, self).__init__()
    self.className = self.__class__.__name__
    self.myName = "%s_%03d" % (self.className, self.noThread[0])
    self.noThread[0] += 1
    self.toStop = False
    self.nbReceive = 0
    self.nbSend = 0
    self._pipeFifo = None
    self._wait = 1. #.5 # 1/2 sec loop of run
    self.verbose = verbose

  def getStrTime(self):
    return time.strftime("%Hh:%Mm:%Ss")

  def run(self):
    if self.verbose:
      logger.info("%s.run() of %s on %s" % (self.className, self.myName, self.getBaseNameFifo()))
    mess = """ThreadPipeFifoBase.run: non implemented in base class:
choose ThreadPipeFifoReceiveBlocking or ThreadPipeFifoReceiveNoBlocking"""
    #raise Exception(mess) #prefer test
    ii = 0
    while True:
      if self.verbose: logger.info("%s: loop for test %i second" % (self.className, ii))
      ii += 1
      time.sleep(1)
      if self.toStop: break
    if self.verbose: logger.info("%s: %s stopped" % (self.className, self.getBaseNameFifo()))

  def stop(self, Verbose=False):
    self.toStop = True
    #and remove file fifoname ???

  def getBaseNameFifo(self):
    if self._pipeFifo == None: return None
    return self._pipeFifo.getBaseNameFifo()

  def getNameFifo(self):
    if self._pipeFifo == None: return None
    return self._pipeFifo.getNameFifo()

  def getPipeFifo(self):
    return self._pipeFifo

  def __repr__(self):
    res = "%s(myName=%s, pipe=%s, nbSend=%i, nbReceive=%i)" % \
          (self.className, self.myName, self.getBaseNameFifo(), self.nbSend, self.nbReceive)
    return res


#######################################################
class ThreadPipeFifoReceiveNoBlocking(ThreadPipeFifoBase):
  """
  Thread execute 'caller.receiveStrSignal.emit(str)' when receive string from fifo
  NoBlocking valid only python
  """

  def __init__(self, caller, fileFifo):
    super(ThreadPipeFifoReceiveNoBlocking, self).__init__()
    self._caller = caller
    self._pipeFifo = PipeFifo(fileFifo)

  def run(self):
    """NoBlocking valid only python"""
    if self.verbose: logger.info("%s.run() of %s on %s" % (self.className, self.myName, self.getBaseNameFifo()))
    while True:
      if self.toStop:
        if verbose: logger.info("%s toStop %s" % (self.myName, self.toStop))
        break
      try:
        message = self._pipeFifo.readNoBlocking()
      except: # case named pipe file deleted by other
        self.stop()
        if self.verbose: logger.warning("ThreadPipeFifoReceiveNoBlocking stop on pipe deleted, or else (%s)" % self.myName)
        continue
      if message != "":
        self.nbReceive += 1
        if self._caller != None:
          if self.verbose:
            logger.info("receiveSignal.emit('%s')" % message)
          self._caller.receiveSignal.emit(message)
        else:
          if self.verbose: logger.warning("no caller, receive: '%s' (%s)" % (message, self.myName))
      time.sleep(self._wait)
    if self.verbose: logger.info('end %s receive run()' % self.myName)


#######################################################
class ThreadPipeFifoReceiveBlocking(ThreadPipeFifoBase):
  """
  Thread execute 'caller.receiveStrSignal.emit(str)' when receive string from fifo
  Blocking valid for python and ROOT (write needs pipe read continously)
  """

  def __init__(self, caller, fileFifo):
    super(ThreadPipeFifoReceiveBlocking, self).__init__()
    self._caller = caller
    self._pipeFifo = PipeFifo(fileFifo)
    self._endPipeMessage = "__endPipe__\n"

  def run(self):
    """blocking valid only python and ROOT"""
    if self.verbose:
      logger.info("%s.run() of %s on %s" % (self.className, self.myName, self.getBaseNameFifo()))
      logger.info("\n!!!! line as readBlocking for ROOT !!!!\n")

    #ii = 0
    while True:
      if self.toStop:
        if verbose: logger.info("%s toStop %s" % (self.myName, self.toStop))
        break
      try:
        message = self._pipeFifo.readlineBlocking()
        #message = "print('coucou %i') #TODO test remove that\n" %ii
        #ii += 1
        #time.sleep(1)
      except: # case named pipe file deleted by other
        self.stop()
        logger.warning("ThreadPipeFifoReceiveBlocking stop on pipe deleted, or else (%s)" % self.myName)
        continue
      if message != "":
        if message == self._endPipeMessage:
          if self.verbose: logger.warning("stop on pipe __endPipe__ (%s)" % self.myName)
          return #as stop
        self.nbReceive += 1
        if self._caller != None:
          if self.verbose: logger.info("receiveSignal.emit('%s')" % message[:-1])
          self._caller.receiveSignal.emit(message)
          #self._caller.receiveMessage(message)
        else:
          if self.verbose: logger.warning("no caller, receive: '%s' (%s)" % (message, self.myName))
      # no wait because blocking read
      # time.sleep(self._wait)
    if self.verbose: logger.info('end %s receive run()' % self.myName)

  def stop(self, Verbose=False):
    self.toStop = True
    aFile = self.getNameFifo()
    try:
      pipeout = os.open(aFile, os.O_WRONLY | os.O_NONBLOCK)   # no blocking if not open "r" at other pipe side
    except OSError as exc: # seems to be when pipe is not open read
      return #nothing to do...
    if Verbose: logger.info("stop %s" % self.getNameFifo())
    os_write(pipeout, self._endPipeMessage)
    os.close(pipeout)


#######################################################
class ThreadPipeFifoSend(ThreadPipeFifoBase):
  """
  Thread execute 'caller.sendStrSignal.emit(str)'  when send string to fifo
  seems immutable objects are automatically threadsafe
  so mutex protecting on string self._message in aThread.setMessage() is useless
  """

  def __init__(self, caller, fileFifo):
    super(ThreadPipeFifoSend, self).__init__()
    self._caller = caller
    self._pipeFifo = PipeFifo(fileFifo)
    self._defaultMessage = "Hello for test %s" % self.className
    self._wait = 1.      # scrute for end thread if PipeFifo deleted, by example

  def sendMessage(self, message):
    if self.is_alive():
      seldMessage = "%s at %s" % (message, self.getStrTime()) #TODO self._message
      ok = self._pipeFifo.writeNoBlocking(message)
      if ok :
        self.nbSend += 1
        if self._caller != None:
          self._caller.endStrSignal.emit(message)
        else:
          if self.verbose: logger.warning("no caller, send:    '%s' (%s)" % (message, self.myName))

    else:
      logger.warning("thread not alive: NO send: '%s' (%s)" % (message, self.myName))

  def run_nonblocking(self):
    """for test"""
    while True:
      if self.toStop:
        if self.verbose: logger.info("%s toStop %s" % (self.myName, self.toStop))
        break
      ok = self._pipeFifo.isPipeFifo() # scrute for end thread if PipeFifo deleted, by example
      if not ok:
        self.stop()
        if self.verbose: logger.warning("ThreadPipeFifoSend stop on pipe deleted, or else (%s)" % self.myName)
        break # stop immediately
      else:
        if self.verbose: logger.warning("pipe ok and ready (%s)" % self.myName)
      time.sleep(self._wait)
    if self.verbose: logger.info('end %s send run_nonblock()' % self.myName)

  def run_nonblocking_test(self):
    """for test to send periodically self._defaultMessage for test"""
    while True:
      if self.toStop:
        if self.verbose: logger.info("%s toStop %s" % (self.myName, self.toStop))
        break
      message = "%s at %s" % (self.message, self.getStrTime()) #TODO self._message
      ok = self._pipeFifo.writeNoBlocking(message)
      if ok :
        self.nbSend += 1
        if self._caller != None:
          self._caller.endStrSignal.emit(message)
        else:
          if self.verbose: logger.warning("no caller, send:    '%s' (%s)" % (message, self.myName))
      time.sleep(self._wait)
    if self.verbose: logger.info('end %s send run_nonblock()' % self.myName)
