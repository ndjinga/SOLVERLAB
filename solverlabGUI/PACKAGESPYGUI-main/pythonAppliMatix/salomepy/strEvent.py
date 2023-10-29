#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


def strEvent(event):
    """catch readable explicit type of event"""
    mapEvent={}
    anInt = event.type()
    #for i in dir(event): print "event", i, getattr(event, i).__class__.__name__, getattr(event, i)
    for i in dir(event):
      anAttr = getattr(event, i)
      if anAttr.__class__.__name__ == "Type":
        #print "event", i, anAttr.__class__.__name__, anAttr
        mapEvent[int(anAttr)] = str(i)
    try:
      res = mapEvent[int(event.type())]
    except:
      res = "unknown event '%i'" % anInt
    #print "strEvent ****** ", event.type(), res
    return res
