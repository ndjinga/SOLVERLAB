#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see SOLVERLABGUI/LICENSE file
# %% LICENSE_END


__all__ = [
    "moduleID",
    "objectID",
    "unknownID",
    "moduleName",
    "modulePixmap",
    "verbose",
    "getORB",
    "getNS",
    "getLCC",
    "getStudyManager",
    "getEngine",
    "getEngineIOR",
    "findOrCreateComponent",
    "getObjectID",
    ]

from omniORB import CORBA
from SALOME_NamingServicePy import SALOME_NamingServicePy_i
from LifeCycleCORBA import LifeCycleCORBA
import SALOMEDS
import SALOMEDS_Attributes_idl
#import SOLVERLABGUI_ORB

_NAME =  "SOLVERLAB"
_NAMEGUI =  "SOLVERLABGUI"
###
# Set/Get verbose level
###
# __verbose__ = int( os.getenv( 'SALOME_VERBOSE', 0 ) ) # for salome trace
# __verbose__ = int( os.getenv( 'SALOME_VERBOSE_MICROGEN', 0 ) ) # for salome MICROGEN trace
__verbose__ = 1  # have to set 0 in production

def log_print(*args):
  if __verbose__ != 0:
    if len(args) == 0:
      return
    if len(args) == 1:
      print(args[0])
    else:
      print(args)

###
# Get module's ID
###
def moduleID():
    MODULE_ID = 71000
    log_print(_NAMEGUI + " MODULE_ID %s" % MODULE_ID)
    return MODULE_ID

###
# Get module object's ID
###
def objectID():
    OBJECT_ID = moduleID() + 10
    return OBJECT_ID

###
# Get unknown ID
###
def unknownID():
    FOREIGN_ID = -1
    return FOREIGN_ID

###
# Get module's name
###
def moduleName():
    return _NAMEGUI

###
# Get module's pixmap name
###
def modulePixmap():
    return _NAMEGUI + "_small.png"

###
# Get ORB reference
###
__orb__ = None
def getORB():
    global __orb__
    if __orb__ is None:
        __orb__ = CORBA.ORB_init( [''], CORBA.ORB_ID )
        pass
    return __orb__

###
# Get naming service instance
###
__naming_service__ = None
def getNS():
    global __naming_service__
    if __naming_service__ is None:
        __naming_service__ = SALOME_NamingServicePy_i( getORB() )
        pass
    return __naming_service__

##
# Get life cycle CORBA instance
##
__lcc__ = None
def getLCC():
    global __lcc__
    if __lcc__ is None:
        __lcc__ = LifeCycleCORBA( getORB() )
        pass
    return __lcc__

##
# Get study manager
###
__study_manager__ = None
def getStudyManager():
    global __study_manager__
    if __study_manager__ is None:
        obj = getNS().Resolve( '/myStudyManager' )
        __study_manager__ = obj._narrow( SALOMEDS.StudyManager )
        pass
    return __study_manager__

###
# Get engine
###
__engine__ = None
def getEngine():
    global __engine__
    if __engine__ is None:
        __engine__ = getLCC().FindOrLoadComponent( "FactoryServerPy", moduleName() )
        pass
    return __engine__

###
# Get engine IOR
###
def getEngineIOR():
    IOR = ""
    if getORB() and getEngine():
        IOR = getORB().object_to_string( getEngine() )
        pass
    return IOR

###
# Find or create module component object in a study
###
def findOrCreateComponent( study ):
    father = study.FindComponent( moduleName() )
    if father is None:
        builder = study.NewBuilder()
        father = builder.NewComponent( moduleName() )
        attr = builder.FindOrCreateAttribute( father, "AttributeName" )
        attr.SetValue( moduleName() )
        attr = builder.FindOrCreateAttribute( father, "AttributePixMap" )
        attr.SetPixMap( modulePixmap() )
        attr = builder.FindOrCreateAttribute( father, "AttributeLocalID" )
        attr.SetValue( moduleID() )
        try:
            builder.DefineComponentInstance( father, getEngine() )
            pass
        except:
            pass
        pass
    return father

###
# Get object's ID
###
def getObjectID( study, entry ):
    ID = unknownID()
    if study and entry:
        sobj = study.FindObjectID( entry )
        if sobj is not None:
            test, anAttr = sobj.FindAttribute( "AttributeLocalID" )
            if test: ID = anAttr._narrow( SALOMEDS.AttributeLocalID ).Value()
            pass
        pass
    return ID
    
