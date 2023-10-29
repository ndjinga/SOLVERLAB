#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
simple tagging as '<color>' for simple coloring log messages on terminal(s)
window or unix or ios using backend colorama.
Using '<color>' because EZ human readable,
So '<color>' are not supposed existing in log message.
"{}".format() is not choosen because "{}" are present
in log messages of contents of python dict (as JSON) etc.

| Usage:
| >> import solverlabpy.coloringSvl as COLS
| 
| Example:
| >> log("this is in <green>color green<reset>, OK is in blue: <blue>OK?")
"""

import os
import sys
import pprint as PP

_verbose = True
_name = "coloringSvl"

""" 
https://github.com/tartley/colorama
init(wrap=True):

On Windows, colorama works by replacing sys.stdout and sys.stderr 
with proxy objects, which override the .write() method to do their work. 
If this wrapping causes you problems, 
then this can be disabled by passing init(wrap=False). 
The default behaviour is to wrap if autoreset or strip or convert are True.

When wrapping is disabled, colored printing on non-Windows platforms 
will continue to work as normal. 
To do cross-platform colored output, 
you can use Colorama's AnsiToWin32 proxy directly:

| Example:
| >>  import sys
| >>  from colorama import init, AnsiToWin32, Fore
| >>  init(wrap=False)
| >>  stream = AnsiToWin32(sys.stderr).stream
| >>  # Python 2
| >>  print >>stream, Fore.BLUE + 'blue text on stderr'
| >>  # Python 3
| >>  print(Fore.BLUE + 'blue text on stderr', file=stream)
"""

import platform
import xyzpy.colorama as CLRM
from solverlabpy.colorama import Fore as FG
from solverlabpy.colorama import Style as ST
#from solverlabpy.colorama import AnsiToWin32
from solverlabpy.colorama import AnsiToWin32 # debug is os.name == 'nt' ?

verbose = False

CLRM.init(wrap=False) # choose NO wrapping

"""
from colorama:
Available formatting constants are:
Fore: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET.
Back: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET.
Style: DIM, NORMAL, BRIGHT, RESET_ALL

note: DIM is not assumed in win32
"""

# dir(ST)

# order matters for items replaces forward to color
_tags = [
  ("<black>", FG.BLACK),
  ("<red>", FG.RED),
  ("<green>", FG.GREEN),
  ("<yellow>", FG.YELLOW),
  ("<blue>", FG.BLUE),
  ("<magenta>", FG.MAGENTA),
  ("<cyan>", FG.CYAN),
  ("<white>", FG.WHITE),
  ("<bright>", ST.BRIGHT),
  ("<normal>", ST.NORMAL),
  ("<reset>", ST.RESET_ALL),
  ("<info>", FG.MAGENTA),
  ("<header>", FG.BLUE),
  ("<label>", FG.CYAN),
  ("<success>", FG.GREEN),
  ("<warning>", FG.RED),
  ("<error>", FG.RED + ST.BRIGHT),
  ("<critical>", FG.RED + ST.BRIGHT),
]

# _tagsNone = ( (i, "") for i,j in _tags ) # to clean tags when log not tty
# reversed order matters for item replaces backward to no color
_tagsNone = list( reversed( [(i, "") for i, j in _tags] ) )

# more non empty colored smart tags
_tags = _tags + [
  ("<OK>", FG.GREEN + ST.BRIGHT + "OK" + ST.RESET_ALL),
  ("<KO>", FG.RED + ST.BRIGHT + "KO" + ST.RESET_ALL),
]

# more non empty colored smart tags reversed order
_tagsNone = [
  ("<OK>", "OK"),
  ("<KO>", "KO"),
] + _tagsNone

# import pprint as PP
# print("_tags %s" % PP.pformat(_tags))
# print("_tagsNone %s" % PP.pformat(_tagsNone))

def indent(msg, nb, car=" "):
  """indent nb car (spaces) multi lines message except first one"""
  s = msg.split("\n")
  res = ("\n"+car*nb).join(s)
  return res

def log(msg):
  """elementary log stdout for debug if _verbose"""
  prefix = "%s.log: " % _name
  nb = len(prefix)
  if _verbose: 
    ini = prefix + indent(msg, nb)
    res = toColor(ini)
    if res != ini:
      res = res + toColor("<reset>")
    print(res)
  
class ColoringStream(object):
  """
  write my stream class
  only write and flush are used for the streaming
  https://docs.python.org/2/library/logging.handlers.html
  https://stackoverflow.com/questions/31999627/storing-logger-messages-in-a-string
  """
  def __init__(self):
    self.logs = ''

  def write(self, astr):
    # log("UnittestStream.write('%s')" % astr)
    self.logs += astr

  def flush(self):
    pass

  def __str__(self):
    return self.logs

def toColor(msg):
  """
  automatically clean the message of color tags '<red> ... 
  if the terminal output stdout is redirected by user
  if not, replace tags with ansi color codes
  """
  if verbose:
    import debogpy.debug as DBG # avoid cross import
    DBG.write("toColor isatty %s %s" % ('isatty' in dir(sys.stdout), sys.stdout.isatty()), msg[0:40])
  if not ('isatty' in dir(sys.stdout) and sys.stdout.isatty()):
    # clean the message color (if the terminal is redirected by user)
    return replace(msg, _tagsNone)
  elif platform.system() == "Windows":
      # clean the message color (is there color in Anaconda prompt or cmd.exe ???)
      return replace(msg, _tagsNone)
  else:
    # https://en.wikipedia.org/wiki/ANSI_escape_code#Escape_sequences
    # rc + CSI 0 K as clear from cursor to the end of the line
    s = msg.replace("<RC>", "\r\033[0;K")
    return replace(s, _tags)
    
def cleanColors(msg):
  """clean the message of color tags '<red> ... """
  return replace(msg, _tagsNone)
  
def toColor_AnsiToWin32(msg):
  """for test debug no wrapping"""
  if not ('isatty' in dir(sys.stdout) and sys.stdout.isatty()):
    # clean the message color if the terminal is redirected by user
    return replace(msg, _tagsNone)
  else:
    msgAnsi = replace(msg, _tags)
    streamOut = ColoringStream()
    atw = AnsiToWin32(streamOut, convert=True)
    streamIn = atw.stream
    print("should_wrap", atw.should_wrap(), atw.convert, atw.strip, atw.autoreset)
    streamIn.write(msgAnsi)
    #AnsiToWin32(streamOut).write_and_convert(msgAnsi)
    # print("streamOut %s" % str(streamOut))
    return str(streamOut)

def replace(msg, tags):
  s = msg
  for r in tags:
    s = s.replace(*r)
  return s
  
if __name__ == "__main__":  
  #log(FG.BLUE + 'blue text on stdout'+ ST.RESET_ALL)
  log("import <green>colorama at <blue>%s" % CLRM.__file__)
  log("import <green>colorama<reset> in <blue>%s: <OK>" % __file__)
  log("import <green>colorama<reset> in <blue>%s: <KO>" % __file__)
  log("import <green>colorama in <blue>%s" % __file__)
  log("set <green>green and not reset...")
  log("...and here is not green because appended reset at end of message")
  log("dir(FG):\n<blue>%s ... <OK> or <KO> ??" % dir(FG))
  log("dir(ST):\n<blue>%s ... <OK> or <KO> ??" % dir(ST))

  
