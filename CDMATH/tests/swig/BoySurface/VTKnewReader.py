#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : Anthony GEAY (CEA/DEN/DM2S/STMF)

from MEDLoader import *

class PVDReader:
    @classmethod
    def New(cls,fileName):
        """ Static constructor. """
        return PVDReader(fileName)
        pass

    def __init__(self,fileName):
        self._fileName=fileName
        pass

    def loadTopInfo(self):
        fd=open(self._fileName,"r")
        return self.__parseXML(fd)

    def __parseXML(self,fd):
        import xml.sax
        class PVD_SAX_Reader(xml.sax.ContentHandler):
            def __init__(self):
                self._tsteps=[]
                pass
            def startElement(self,name,attrs):
                if name=="VTKFile":
                    if attrs["type"]!="Collection":
                        raise Exception("Mismatch between reader (PVD) type and file content !")
                    return
                if name=="DataSet":
                    self._tsteps.append((float(attrs["timestep"]),str(attrs["file"])))
                    return
                pass
            pass
        rd=PVD_SAX_Reader()
        parser=xml.sax.make_parser()
        parser.setContentHandler(rd)
        parser.parse(fd)
        return rd
    pass

class PVTUReader:
    @classmethod
    def New(cls,fileName):
        """ Static constructor. """
        return PVTUReader(fileName)
        pass

    def __init__(self,fileName):
        self._fileName=fileName
        pass

    def loadParaInfo(self):
        fd=open(self._fileName,"r")
        return self.__parseXML(fd)

    def __parseXML(self,fd):
        import xml.sax
        class PVTU_SAX_Reader(xml.sax.ContentHandler):
            def __init__(self):
                self._data_array={2:self.DAPointData,3:self.DACellData}
                self._node_fields=[]
                self._cell_fields=[]
                self._pfiles=[]
                self._tmp=None
                pass
            def DAPointData(self,attrs):
                self._node_fields.append((str(attrs["Name"]),int(attrs["NumberOfComponents"])))
                pass
            def DACellData(self,attrs):
                self._cell_fields.append((str(attrs["Name"]),int(attrs["NumberOfComponents"])))
                pass
            def startElement(self,name,attrs):
                if name=="VTKFile":
                    if attrs["type"]!="PUnstructuredGrid":
                        raise Exception("Mismatch between reader (PVTU) type and file content !")
                    return
                if name=="Piece":
                    self._pfiles.append(str(attrs["Source"]))
                    return
                if name=="PPointData":
                    self._tmp=2
                    return
                if name=="PCellData":
                    self._tmp=3
                    return
                if name=="PDataArray":
                    if self._tmp in self._data_array:
                        self._data_array[self._tmp](attrs)
                        pass
                    return
                pass
            pass
        rd=PVTU_SAX_Reader()
        parser=xml.sax.make_parser()
        parser.setContentHandler(rd)
        parser.parse(fd)
        return rd
    pass

class VTURawReader:
    """ Converting a VTU file in raw mode into the MED format.
    Warning: VTU file must be write in "Appended" mode, and in "Binary" format.
    """
    VTKTypes_2_MC=[-1,0,-1,1,33,NORM_TRI3 ,-1,5,-1,4,14,-1,NORM_HEXA8,16,15,-1,22,-1,-1,-1,-1,2,6,8,20,30,25,23,9,27,-1,-1,-1,-1,7,-1,-1,-1,-1,-1,-1,-1,31]

    class NormalException(Exception):
        def __init__(self,lineNb):
            Exception.__init__(self)
            self._line_nb=lineNb
        def getLineNb(self):
            return self._line_nb
        pass

    class NotRawVTUException(Exception):
        pass

    def loadInMEDFileDS(self):
        import numpy as np
        fd=open(self._fileName,"r")
        ref,rd=self.__parseXML(fd)
        #
        ret=MEDFileData()
        ms=MEDFileMeshes() ; ret.setMeshes(ms)
        fs=MEDFileFields() ; ret.setFields(fs)
        #
        types=np.memmap(fd,dtype=rd._type_types,mode='r',offset=ref+rd._off_types,shape=(rd._nb_cells,))
        types=self.__swapIfNecessary(rd._bo,types)
        # mesh dimension detection
        types2=types.copy() ; types2.sort() ; 
        types2 =np.unique(types2)
        # Get first valid type
        meshDim = -1
        for i, typ in enumerate(types2):
            if self.VTKTypes_2_MC[typ] != -1:
                meshDim=MEDCouplingMesh.GetDimensionOfGeometricType(self.VTKTypes_2_MC[typ])
        if meshDim == -1:
            raise Exception("Could not find a valid cell type in the mesh !")
        if i > 0:
            print("WARNING: some invalid/incompatible cell types were detected - trying to have them as polygons ...")
        for typ in types2:
            if self.VTKTypes_2_MC[typ] != -1:
                md=MEDCouplingMesh.GetDimensionOfGeometricType(self.VTKTypes_2_MC[typ])
                if md!=meshDim:
                    raise Exception("MultiLevel umeshes not managed yet !")
            else:
                print("WARNING: invalid/incompatible cell type detected: VTK type number: %d" % typ)
        m=MEDCouplingUMesh("mesh",meshDim)
        # coordinates
        coo=np.memmap(fd,dtype=rd._type_coords,mode='r',offset=ref+rd._off_coords,shape=(rd._nb_nodes*rd._space_dim,))
        coo=self.__swapIfNecessary(rd._bo,coo) ; coo=DataArrayDouble(np.array(coo,dtype='float64')) ; coo.rearrange(rd._space_dim)
        m.setCoords(coo)
        # connectivity
        offsets=np.memmap(fd,dtype=rd._type_off,mode='r',offset=ref+rd._off_off,shape=(rd._nb_cells,))
        offsets=self.__swapIfNecessary(rd._bo,offsets) ; connLgth=offsets[-1] ; offsets2=DataArrayIdType(rd._nb_cells+1) ; offsets2.setIJ(0,0,0)
        offsets2[1:]=DataArrayIdType([int(o) for o in offsets])
        offsets3=offsets2.deltaShiftIndex() ; offsets2=offsets3.deepCopy() ; offsets3+=1 ; offsets3.computeOffsetsFull()
        offsets=offsets3
        tmp1=DataArrayIdType(len(offsets2),2) ; tmp1[:,0]=1 ; tmp1[:,1]=offsets2 ; tmp1.rearrange(1) ; tmp1.computeOffsetsFull()
        tmp1=DataArrayIdType.Range(1,2*len(offsets2),2).buildExplicitArrByRanges(tmp1)
        conn=np.memmap(fd,dtype=rd._type_conn,mode='r',offset=ref+rd._off_conn,shape=(connLgth,))
        conn=self.__swapIfNecessary(rd._bo,conn)
        types=np.array(types,dtype='int32') ; types=DataArrayIdType(types) ; 
        types.transformWithIndArr(self.VTKTypes_2_MC)
        conn2=DataArrayIdType(offsets.back())
        conn2[offsets[0:-1]]=types
        conn2[tmp1]=DataArrayIdType([int(c) for c in conn])
        m.setConnectivity(conn2,offsets,True)
        m.checkConsistencyLight() ; mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; ms.pushMesh(mm)
        # Fields on nodes and on cells
        for spatialDisc,nbEnt,fields in [(ON_NODES,rd._nb_nodes,rd._node_fields),(ON_CELLS,rd._nb_cells,rd._cell_fields)]: 
            for name,typ,nbCompo,off in fields:
                ff=MEDFileFieldMultiTS()
                f=MEDCouplingFieldDouble(spatialDisc,ONE_TIME)
                f.setName(name) ; f.setMesh(m)
                vals=np.memmap(fd,dtype=typ,mode='r',offset=ref+off,shape=(nbEnt*nbCompo))
                vals=self.__swapIfNecessary(rd._bo,vals)
                arr=DataArrayDouble(np.array(vals,dtype='float64')) ; arr.rearrange(nbCompo)
                f.setArray(arr) ; f.checkConsistencyLight()
                f.setTime(self._time[0],self._time[1],0)
                ff.appendFieldNoProfileSBT(f)
                fs.pushField(ff)
                pass
            pass
        return ret

    def __parseXML(self,fd):
        import xml.sax
        class VTU_SAX_Reader(xml.sax.ContentHandler):
            def __init__(self):
                self._loc=None
                self._data_array={0:self.DAPoints,1:self.DACells,2:self.DAPointData,3:self.DACellData}
                self._node_fields=[]
                self._cell_fields=[]
                pass
            def setLocator(self,loc):
                self._loc=loc
            def DAPoints(self,attrs):
                self._space_dim=int(attrs["NumberOfComponents"])
                self._type_coords=str(attrs["type"]).lower()
                self._off_coords=int(attrs["offset"])
                pass
            def DACells(self,attrs):
                if attrs["Name"]=="connectivity":
                    self._type_conn=str(attrs["type"]).lower()
                    self._off_conn=int(attrs["offset"])
                    pass
                if attrs["Name"]=="offsets":
                    self._type_off=str(attrs["type"]).lower()
                    #self._off_off=int(attrs.get("offset", 0))
                    self._off_off=int(attrs["offset"])
                    pass
                if attrs["Name"]=="types":
                    self._type_types=str(attrs["type"]).lower()
                    #self._off_types=int(attrs.get("offset", 0))
                    self._off_types=int(attrs["offset"])
                    pass
                pass
            def DAPointData(self,attrs):
                numCompo = int(attrs.get("NumberOfComponents", 1))
                offset = int(attrs.get("offset", 0))
                self._node_fields.append((str(attrs["Name"]),str(attrs["type"]).lower(),numCompo,offset))
                pass
            def DACellData(self,attrs):
                self._cell_fields.append((str(attrs["Name"]),str(attrs["type"]).lower(),int(attrs["NumberOfComponents"]),int(attrs["offset"])))
                pass
            def startElement(self,name,attrs):
                if name=="VTKFile":
                    if attrs["type"]!="UnstructuredGrid":
                        raise Exception("Mismatch between reader (VTU) type and file content !")
                    self._bo=bool(["LittleEndian","BigEndian"].index(attrs["byte_order"]))
                    pass
                if name=="Piece":
                    self._nb_cells=int(attrs["NumberOfCells"])
                    self._nb_nodes=int(attrs["NumberOfPoints"])
                    return
                if name=="Points":
                    self._tmp=0
                    return
                if name=="Cells":
                    self._tmp=1
                    return
                if name=="PointData":
                    self._tmp=2
                    return
                if name=="CellData":
                    self._tmp=3
                    return
                if name=="DataArray":
                    self._data_array[self._tmp](attrs)
                    return
                if name=="AppendedData":
                    if str(attrs["encoding"])=="raw":
                        raise VTURawReader.NormalException(self._loc.getLineNumber())
                    else:
                        print(attrs["encoding"])
                        raise VTURawReader.NotRawVTUException("The file is not a raw VTU ! Change reader !")
                pass
            pass
        rd=VTU_SAX_Reader()
        parser=xml.sax.make_parser()
        parser.setContentHandler(rd)
        locator=xml.sax.expatreader.ExpatLocator(parser)
        rd.setLocator(locator)
        isOK=False
        try:
            parser.parse(fd)
        except self.NormalException as e:
            isOK=True
            fd.seek(0)
            for i in range(e.getLineNb()): 
              l = fd.readline()
            ref=fd.tell()+12
            pass
        if not isOK:
            raise Exception("Error in VTURawReader : not a raw format ?")
        return ref,rd

    @classmethod
    def New(cls,fileName,tim=(0.,0)):
        """ Static constructor. """
        return VTURawReader(fileName,tim)
        pass

    def __init__(self,fileName,tim=(0.,0)):
        msg="The time specified in constructor as 2nd arg should be a tuple containing 2 values 1 float and 1 int !"
        if not isinstance(tim, tuple):
            raise Exception(msg)
        if len(tim)!=2:
            raise Exception(msg)
        if not isinstance(tim[0], float) or not isinstance(tim[1], int):
            raise Exception(msg)
        self._fileName=fileName
        self._time=tim
        pass

    def __swapIfNecessary(self,b,arr):
        if b:
            ret=arr.copy()
            ret.byteswap(True)
            return ret
        else:
            return arr
        pass
    pass
