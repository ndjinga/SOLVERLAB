#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END


import SOLVERLABGUI_ORB__POA
import SALOME_ComponentPy
import SALOME_DriverPy

from SOLVERLABGUI_utils import *

class SOLVERLABGUI(SOLVERLABGUI_ORB__POA.SOLVERLABGUI_Gen,
              SALOME_ComponentPy.SALOME_ComponentPy_i,
              SALOME_DriverPy.SALOME_DriverPy_i):
    """
    Construct an instance of SOLVERLABGUI module engine.
    The class SOLVERLABGUI implements CORBA interface SOLVERLABGUI_Gen (see SOLVERLABGUI_Gen.idl).
    It is inherited from the classes SALOME_ComponentPy_i (implementation of
    Engines::EngineComponent CORBA interface - SALOME component) and SALOME_DriverPy_i
    (implementation of SALOMEDS::Driver CORBA interface - SALOME module's engine).
    """
    def __init__ ( self, orb, poa, contID, containerName, instanceName, 
                   interfaceName ):
        super().__init__(orb, poa,
                    contID, containerName, instanceName, interfaceName, 0)
        super().__init__(interfaceName)
        self._naming_service = SALOME_ComponentPy.SALOME_NamingServicePy_i( self._orb )
 
    """
    Generate banner.
    """
    def makeBanner( self, name ):
        banner = "Hello %s!" % name
        return banner

    """
    Create object.
    """
    def createObject( self, study, name ):
        builder = study.NewBuilder()
        father  = findOrCreateComponent( study )
        object  = builder.NewObject( father )
        attr    = builder.FindOrCreateAttribute( object, "AttributeName" )
        attr.SetValue( name )
        attr    = builder.FindOrCreateAttribute( object, "AttributeLocalID" )
        attr.SetValue( objectID() )

    """
    Dump module data to the Python script.
    """
    def DumpPython( self, study, isPublished ):
        abuffer = []
        abuffer.append( "def RebuildData( theStudy ):" )
        names = []
        father = study.FindComponent( moduleName() )
        if father:
            iter = study.NewChildIterator( father )
            while iter.More():
                name = iter.Value().GetName()
                if name: names.append( name )
                iter.Next()
        if names:
            abuffer += [ "  from batchmode_salome import lcc" ]
            abuffer += [ "  import SOLVERLABGUI_ORB" ]
            abuffer += [ "  " ]
            abuffer += [ "  myCompo = lcc.FindOrLoadComponent( 'FactoryServerPy', '%s' )" % moduleName() ]
            abuffer += [ "  " ]
            abuffer += [ "  myCompo.createObject( theStudy, '%s' )" % name for name in names ]
        abuffer += [ "  " ]
        abuffer.append( "  pass" )
        abuffer.append( "\0" )
        return ("\n".join( abuffer ), 1)
