#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2010-2018  CEA/DEN
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

"""
utilities for solverlabGUI
general useful simple methods
all-in-one import solverlabpy.utilsSvl as UTS

| Usage:
| >> import solverlabpy.utilsSvl as UTS
| >> UTS.ensure_path_exists(path)
"""

import os
import shutil
import errno
import stat
import time

import re
import tempfile
import subprocess as SP

import returncodepy.returnCode as RCO
import debogpy.debug as DBG # Easy print stderr (for DEBUG only)


##############################################################################
# file system utilities
##############################################################################
def ensure_path_exists(path):
    """Create a path if not existing
    
    :param path: (str) The path.
    """
    # DBG.write("ensure_path_exists", path, True)
    if not os.path.exists(path):
        os.makedirs(path)
        
def ensure_file_exists(aFile, aDefaultFile):
    """
    Create a file if not existing,
    copying from default file
    
    :param aFilepath: (str) The file to ensure existence
    :param aDefaultFile: (str) The default file to copy if not existing
    """
    isfile = os.path.isfile(aFile)
    if isfile: return True
    try:
      DBG.write("ensure_file_exists %s" % isfile, aDefaultFile + " -->\n" + aFile)
      shutil.copy2(aDefaultFile, aFile)
      return True
    except:
      return False

        
def replace_in_file(file_in, str_in, str_out):
    """
    Replace <str_in> by <str_out> in file <file_in>.
    save a file old version as file_in + '_old'

    :param file_in: (str) The file name
    :param str_in: (str) The string to search
    :param str_out: (str) The string to replace.    
    """
    with open(file_in, "r") as f: 
      contents = f.read()
    shutil.move(file_in, file_in + "_old")
    with open(file_in, "w") as f: 
      f.write(contents.replace(str_in, str_out))


def get_solverlabGUI_version(config):
    # CVW TODO
    return version

def get_tmp_filename(config, name):
    # CVW TODO
    return os.path.join(config.VARS.tmp_root, name)

##############################################################################
# logger utilities
##############################################################################
def formatTuples(tuples):
    """
    Format 'label = value' the tuples in a tabulated way.
    
    :param tuples: (list) The list of tuples to format
    :return: (str) The tabulated text. (as mutiples lines)
    """
    # find the maximum length of the first value of the tuples
    smax = max([len(l[0]) for l in tuples])
    # Print each item of tuples with good indentation
    msg = ""
    for i in tuples:
        sp = " " * (smax - len(i[0]))
        msg += sp + "%s = %s\n" % (i[0], i[1]) # tuples, may be longer
    if len(tuples) > 1: msg += "\n" # skip one line for long list
    return msg
    
def formatValue(label, value, suffix=""):
    """
    format 'label = value' with the info color
    
    :param label: (int) the label to print.
    :param value: (str) the value to print.
    :param suffix: (str) the optionnal suffix to add at the end.
    """
    msg = "  %s = %s %s" % (label, value, suffix)
    return msg


##############################################################################
# color utilities, for convenience    
##############################################################################
_colors = "BLACK RED GREEN YELLOW BLUE MAGENTA CYAN WHITE".lower().split(" ")
    
def black(msg):
    return "<black>"+msg+"<reset>"

def red(msg):
    return "<red>"+msg+"<reset>"

def green(msg):
    return "<green>"+msg+"<reset>"

def yellow(msg):
    return "<yellow>"+msg+"<reset>"

def blue(msg):
    return "<blue>"+msg+"<reset>"

def magenta(msg):
    return "<magenta>"+msg+"<reset>"

def cyan(msg):
    return "<cyan>"+msg+"<reset>"

def white(msg):
    return "<white>"+msg+"<reset>"

def normal(msg):
    return "<normal>"+msg+"<reset>"

def reset(msg):
    return "<reset>"+msg

def info(msg):
    return "<info>"+msg+"<reset>"

def header(msg):
    return "<info>"+msg+"<reset>"

def label(msg):
    return "<label>"+msg+"<reset>"

def success(msg):
    return "<success>"+msg+"<reset>"

def warning(msg):
    return "<warning>"+msg+"<reset>"

def error(msg):
    return "<error>"+msg+"<reset>"

def critical(msg):
    return "<critical>"+msg+"<reset>"
  
def tabColor(*args):
    """
    return tabulated colored string from args, 
    assume true length of color tags as <OK> <info> etc. to correct alignment
    when tags are interpreted as (no-length-spacing) color
    for colorama use or else
    """
    # DBG.write("tabulate args %s" % type(args), args, True)
    i = 0
    imax = len(args)
    res = ""
    for ii in range(30): # no more 30 items
      idx, aStr = args[i:i+2]
      if idx > 0: 
        res += aStr + addSpaces(idx, aStr) # left aligment
      else:
        res += addSpaces(idx, aStr) + aStr # right aligment
      i += 2
      if i >= imax:
        return res
    # something wrong
    #Raise Exception("tabColor problem on %s" % args)
      
def addSpaces(idx, aStr):
    import solverlabpy.coloringSvl as COLS
    cleaned = COLS.cleanColors(aStr)
    # DBG.write("cleaned", "'%s' ->\n'%s'" % (aStr, cleaned), True)
    lg = abs(idx) - len(cleaned)
    if lg <= 0: 
      return ""
    else:
      return " "*lg

##############################################################################
# list and dict utilities
##############################################################################
def deepcopy_list(input_list):
    """Do a deep copy of a list
    
    :param input_list: (list) The list to copy
    :return: (list) The copy of the list
    """
    res = []
    for elem in input_list:
        res.append(elem)
    return res

def remove_item_from_list(input_list, item):
    """Remove all occurences of item from input_list
    
    :param input_list: (list) The list to modify
    :return: (list) The without any item
    """
    res = []
    for elem in input_list:
        if elem == item:
            continue
        res.append(elem)
    return res

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
    


##############################################################################
# subprocess utilities, with logger functionalities (trace etc.)
##############################################################################
    
def Popen(command, shell=True, cwd=None, env=None, stdout=SP.PIPE, stderr=SP.PIPE, logger=None):
  """
  make subprocess.Popen(cmd), with 
  call logger.trace and logger.error if problem as returncode != 0 
  """
  if True: #try:  
    proc = SP.Popen(command, shell=shell, cwd=cwd, env=env, stdout=stdout, stderr=SP.STDOUT)
    res_out, res_err = proc.communicate() # res_err = None as stderr=SP.STDOUT
    rc = proc.returncode

    res_out = res_out.decode("utf-8") # python3 stdout is b'...'
    DBG.write("Popen logger returncode", (rc, res_out))
    
    if rc == 0:
      if logger is not None:
        logger.info("<OK> launch command rc=%s cwd=<info>%s<reset>:\n%s" % (rc, cwd, command))
        logger.info("<OK> result command stdout&stderr:\n%s" % res_out)
      return RCO.ReturnCode("OK", "Popen command done", value=res_out)
    else:
      if logger is not None:
        logger.warning("<KO> launch command rc=%s cwd=<info>%s<reset>:\n%s" % (rc, cwd, command))
        logger.warning("<KO> result command stdout&stderr:\n%s" % res_out)
      return RCO.ReturnCode("KO", "Popen command problem", value=res_out)
  else: #except Exception as e:
    logger.error("<KO> launch command cwd=%s:\n%s" % (cwd, command))
    logger.error("launch command exception:\n%s" % e)
    return RCO.ReturnCode("KO", "Popen command problem")


def sleep(sec):
    time.sleep(sec)