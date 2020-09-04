import pvsimple as pvs
from pvsimple import *
import os
import time
import numpy as np
import MEDReader
from vtk.util import numpy_support as npvtk 
# do I need to kill the pipeline?

def Extract_MED_data_over_line_to_csv_file(inputFileName, outputFileName,
                                           point1, point2,
                                           resolution
                                           ):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = MEDReader.Reader(FileName=[inputFileName])
    PlotOverLine1 = pvs.PlotOverLine(Source="High Resolution Line Source"
                                     )
    PlotOverLine1.Source.Point1 = point1
    PlotOverLine1.Source.Point2 = point2
    PlotOverLine1.Source.Resolution = resolution
    writer = pvs.CreateWriter(outputFileName, PlotOverLine1)
    writer.FieldAssociation = "Cells"  # or "Points" (default in cdmath is finite volumes)
    writer.UpdatePipeline()
    
def Extract_MED_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = MEDReader.Reader(FileName=[inputFileName])
    PlotOverLine1 = pvs.PlotOverLine(Source="High Resolution Line Source"
                                     )
    PlotOverLine1.Source.Point1 = point1
    PlotOverLine1.Source.Point2 = point2
    PlotOverLine1.Source.Resolution = resolution

    vtkarray = PlotOverLine1.GetCellData() # or Slice1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array
    
def Extract_field_data_over_line_to_numpyArray(field, point1, point2, resolution):
    field.writeMED(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+".med"

    result = Extract_MED_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
    return result

def Slice_MED_data_over_line_to_numpyArray(inputFileName, point, normal):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = MEDReader.Reader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    

    vtkarray = Slice1.GetCellData() # or PlotOverLine1.GetCellData() # or Clip1.GetCellData()
    numpy_array = npvtk.vtk_to_numpy(vtkarray)

    return numpy_array

def Slice_MED_data_over_line_to_csv_file(inputFileName,
                                        outputFileName,
                                        point,
                                        normal):
    pvs._DisableFirstRenderCameraReset()
    data_vtu = MEDReader.Reader(FileName=[inputFileName])
    Slice1 = pvs.Slice(SliceType="Plane")
    Slice1.SliceOffsetValues = [0.0]
    Slice1.SliceType.Origin = point
    Slice1.SliceType.Normal = normal
    CellCenters1 = pvs.CellCenters()    
    writer = pvs.CreateWriter(outputFileName, CellCenters1)
    writer.Precision=30
    writer.FieldAssociation = "Cells"  # or "Points" (default in cdmath is finite volumes)
    writer.UpdatePipeline()
 
def Slice_field_data_over_line_to_numpyArray(field,
                                        outputFileName,
                                        point,
                                        normal):
    field.writeMED(field.get_name())
    inputFileName = os.getcwd()+field.get_name()+".med"
 
    result = Slice_MED_data_over_line_to_numpyArray(inputFileName, point1, point2, resolution)

    os.remove(inputFileName)
    return result

def Save_MED_to_picture_file(inputFileName, field_name,POINTS_or_CELLS,
                             outputFileName,
                             ):
    data_med = MEDReader.Reader(FileName=[inputFileName])
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1052, 476]

    # show data in view
    display = Show(data_med, renderView1)

    # reset view to fit data
    renderView1.ResetCamera()

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

    # set scalar coloring
    ColorBy(display, (POINTS_or_CELLS, field_name))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'Initialvariablesforsphericalexplosion'
    resultfieldLUT = GetColorTransferFunction(field_name)

    # get opacity transfer function/opacity map for 'Initialvariablesforsphericalexplosion'
    resultfieldPWF = GetOpacityTransferFunction(field_name)

    SaveScreenshot(outputFileName) 
