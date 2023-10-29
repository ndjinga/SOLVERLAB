#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os

def readSimpleTableFile( afile ):
    """
    read a simple file with float values in columns, lines:
    
    - first line headers beginning with '#'
      number of colums and titlesxy are extracted
    - first line no headers, not beginning with '#'
      titlesxy =['col0','col1',...]
    - separator blanks.
    - multiple consecutive blanks as only one.
    - this method is not designed for big files.
    """
    aDir, aName = os.path.split(afile)
    if not os.path.exists(afile): return None, aName
    with open(afile, "r") as F: s=F.readlines()
    
    #check if first line begins with #, extract titles_xy
    line = ' '.join(s[0].split()).split()
    firsttag = line[0]
    if firsttag[0] == "#":
      if len(firsttag) == 1:
        titles_xy = line[1:]
      else:
        titles_xy = line
        titles_xy[0] = titles_xy[0][1:] #remove # if '#aName0 aName1...'
      lg = len(titles_xy)
      xy = [[] for i in range(lg)]
      idep = 1
    else: #no header, number of colums is from first line
      lg = len(line)
      titles_xy = ["col%i"%i for i in range(lg)] #default names of columns
      xy = [[] for i in range(lg)]
      idep = 0

    for aline in s[idep:]:
      line = ' '.join(aline.split()).split()
      firsttag = line[0]
      if firsttag[0] == "#": continue #skip other lines beginning with #
      for i in range(lg): #read only headers length lg columns
        xy[i].append(float(line[i].replace("D","e"))) #"1.2D+00" for example
    return xy, titles_xy, aName

