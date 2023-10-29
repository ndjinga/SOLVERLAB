#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
Utils for pandas and OSCAR
"""

import os, sys
import pandas as pd
import numpy as np

verbose = False

def read_oscar_csv(aNameFile):
  """
  read a csv file from OSCAR
  return 2 pandas DataFrame(s):
    infos: from header, 2 columns
    contents: 2d array of data
  """
  
  def is_number(s):
    try:
      float(s)
      return True
    except ValueError:
      pass
    return False
    
  def shortLine(line):
    """trunk too long lines ..."""
    if len(line)>70:
      return line[0:70]+" ..."
    return line
  
  def getHeader(aNameFile, delimiter=";", nbLinesMax=20):
    """
    get length of header and header as dict.
    length is numer of lines to skip before DataFrameModel as csv
    check info in header as first lines (no more nbLinesMax)
    with no more than one or two delimiter
    """
    #print firsts line... check info in header as first lines with on or two delimiter
    lines = []
    if verbose: print("file:",aNameFile)
    with open(aNameFile) as myFile:
      for x in range(nbLinesMax):
        line = next(myFile).strip()
        lineSplit = line.split(delimiter)
        #print "isnumeric?", is_number(lineSplit[0]), lineSplit[0]
        if is_number(lineSplit[0]): #first begin array2d oscar
          break
        if lineSplit[0] == "NomRegion":
          #first exception beginning array2d oscar with no number as first column
          lines.append(line)
          break
        lines.append(line)
        
    if verbose: print("begin array2d: '%s'" % shortLine(lines[-1]))
    aDict = getHeaderDictFromLines(lines, delimiter)
    return aDict
    
  def getIndexedTitles(listTitles):
    res = ""
    i = 0
    for title in listTitles:
      res += "%i: '%s' \n" % (i, title)
      i += 1
    return res
    
  def getHeaderDictFromLines(lines, delimiter=";"):
    """
    return a dict from header of oscar csv file
    first column are keys
    second column are values
    append in dict (for use) some keys:
      initialFile,nbHeader,nbHeaderPandas
      nbCols,columnTitles,columnTitlesMore
      nbheader is number of lines of header
      nbHeaderPandas is number of lines of header excluding empty lines
    """
    aDict = {}
    nbLines = len(lines)
    if verbose: print("getHeaderDictFromLines: nblines", nbLines)
    aDict["nbHeader"] = str(nbLines-1)
    nbHeaderPandas = nbLines-1
    aDict["nbCols"] = str(len(lines[-1].split(delimiter)))

    if nbLines == 1: #no header, columnTitles is first line
      len1 = len(lines[-1].split(delimiter))
      titles = getIndexedTitles(lines[-1].split(delimiter))
      aDict['columnTitles'] = titles
      aDict['columnTitlesMore'] = ""
      aDict["nbHeaderPandas"] = str(nbHeaderPandas)
      return aDict

    #there is header
    for line in lines[:-1]:
      if verbose: print("--> '%s'" % shortLine(line))
      lineSplit = line.split(delimiter)
      if line == "":
        #there is empty line that pandas do not count skipping header in read_csv
        nbHeaderPandas += -1
        pass
      else:
        if lineSplit[0] != "":
          aDict[lineSplit[0]] = lineSplit[1]
          if verbose: print("add aDict '%s': '%s'" % (lineSplit[0], lineSplit[1]))

    len1 = len(lines[-1].split(delimiter))
    len2 = len(lines[-2].split(delimiter))
    titles = getIndexedTitles(lines[-1].split(delimiter))
    aDict['columnTitles'] = titles
    if len2 in [len1, len1-1, len1+1]:
      titles = getIndexedTitles(lines[-2].split(delimiter))
      aDict['columnTitlesMore'] = titles
    else:
      aDict['columnTitlesMore'] = ""
    aDict["nbHeaderPandas"] = str(nbHeaderPandas)
    return aDict


  if aNameFile == '': return
  delimiter = ";"
  aDict = getHeader(aNameFile, delimiter=delimiter)
  aDict["initialFile"] = aNameFile
  nbHeader = int(aDict["nbHeader"])
  nbHeaderPandas = int(aDict["nbHeaderPandas"])
  nbCols = int(aDict["nbCols"])
  if verbose:
    print("info: nbHeader %i nbCols %i" % (nbHeader, nbCols))
    i = 0
    for k in sorted(aDict.keys()):
      print("  %i '%s': '%s'" % (i, k, shortLine(aDict[k])))
      i += 1
  
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  #print "nbHeaderPandas", nbHeader, nbHeaderPandas
  df2 = pd.read_csv(aNameFile, delimiter=delimiter, header=nbHeaderPandas, usecols=list(range(nbCols)))
  return (df1, df2)

def getExampleDataFrameEllipse(periods = 1000):
    df = pd.DataFrame({"teta": np.linspace(0, np.pi*2*5, periods)} )
    #vectorized
    r = np.linspace(1, 0.1, periods)
    df["X"] = np.cos(df.teta)*r
    df["Y"] = np.sin(df.teta)*r
    infos = {"suptitle": "from getExampleDataFrameEllipse", 
             "title": "ellipse", 
             "xtitle": "x", 
             "ytitle": "y=f(x)", 
             "periods": periods, 
             }
    df_infos = pd.DataFrame(infos, index=[0]).transpose()
    """
    aPlot = df.plot(x="X", y="Y", label="ellipse", legend=False)
    lim = 1.5
    aPlot.set_xlim(-lim,lim)
    aPlot.set_ylim(-lim,lim)
    """
    return (df_infos, df)

def getExampleDataFrameElementary(periods = 4):
    df = pd.DataFrame({"X": np.linspace(1, periods, periods)} )
    df["Y"] = df.X*2
    infos = {"suptitle": "from getExampleDataFrameElementary", 
             "title": "y=x*2", 
             "xtitle": "x", 
             "ytitle": "y=f(x)", 
             "periods": periods, 
             }
    df_infos = pd.DataFrame(infos, index=[0]).transpose()
    return (df_infos, df)
    
