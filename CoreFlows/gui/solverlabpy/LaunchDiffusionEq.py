import os

import solverlab

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()


def getFormat(computation):
  if computation.file_format == "MED":
    logger.info("Result file in MED format")
    return solverlab.MED
  logger.info("Result file in VTK format")
  return solverlab.VTK


def getFinite(diffusion):
  if diffusion.finite_method == "Element":
    logger.info("Finite element calculation")
    return (True, solverlab.NODES)
  logger.info("Finite volume calculation")
  return (False, solverlab.CELLS)


def getMethod(numerical):
  if numerical.method == "Explicit":
    logger.info("Method: Explicit")
    return solverlab.Explicit
  logger.info("Method: Implicit")
  return solverlab.Implicit


def getLinear(numerical):
  if numerical.linear_solver == "GMRES":
    logger.info("Linear solver: GMRES")
    return solverlab.GMRES
  logger.info("Linear solver: BICGSTAB")
  return solverlab.BICGSTAB


def getPreconditioner(numerical):
  if numerical.preconditioner == "LU":
    logger.info("Preconditioner: LU")
    return solverlab.LU
  if numerical.preconditioner == "ILU":
    logger.info("Preconditioner: ILU")
    return solverlab.ILU
  logger.info("Preconditioner: None")
  return None


def getMaxIter(numerical):
  return int(numerical.max_iteration)


def setBoundary(boundary,medfile,finite_method,mypb): # TODO not verified
  boundaryField = solverlab.Field(medfile, finite_method, boundary.fieldName, boundary.time_iteration,
                                   boundary.time_sub_iteration, boundary.meshLevel)
  if boundary.type == "Neumann":
    mypb.setNeumannBoundaryCondition(boundary.name,boundaryField)
  else:
    mypb.setDirichletBoundaryCondition(boundary.name,boundaryField)


def launchDEQ(model):
  medfile = os.path.splitext(model.GeometryMed.fileMed)[0]  # remove ".med" solverlab doesn't need it
  n = model.Analysis.caseSolverlab.Equation
  diffusion = model.Model[int(n)]

  logger.info("Preparing to launch the simulation")
  mandatory = diffusion.mandatory_values

  finite_bool, finite_method = getFinite(diffusion)
  # myProblem
  mypb = solverlab.DiffusionEquation(diffusion.space_dim, finite_bool, mandatory.solid_density, mandatory.specific_heat, mandatory.thermal_conductivity)

  # Set the mesh and initial data
  f = diffusion.field_option
  mypb.setInitialField(medfile, diffusion.field_name, f.time_iteration, f.time_sub_iteration, f.meshLevel, finite_method)

  # set optional values
  optional = diffusion.optional_values
  if optional.mode_fluid_temp == "Field":
    f = optional.fluid_temp_field
    fluid_field = solverlab.Field(medfile, finite_method, f.field_name, f.time_iteration, f.time_sub_iteration, f.meshLevel)
    mypb.setFluidTemperatureField(fluid_field)
  else:
    f = optional.fluid_temp_scalar
    mypb.setFluidTemperature(f)

  mypb.setHeatTransfertCoeff(optional.heat_transfer_coeff)

  if optional.mode_heat_power == "Field":
    f = optional.heat_power_field
    heat_field = solverlab.Field(medfile, finite_method, f.field_name, f.time_iteration, f.time_sub_iteration, f.meshLevel)
    mypb.setHeatPowerField(heat_field)
  else:
    f = optional.heat_power_scalar
    mypb.setHeatSource(f)

  if finite_bool:
    boundaryNodeGroupNames = mypb.getMesh().getNameOfNodeGroups()
    logger.info("%s Boundary Node Group detected: %s" % (len(boundaryNodeGroupNames), boundaryNodeGroupNames))
  else:
    boundaryFaceGroupNames = mypb.getMesh().getNameOfFaceGroups()
    logger.info("%s Boundary Face Group detected: %s" % (len(boundaryFaceGroupNames), boundaryFaceGroupNames))

  boundary = diffusion.boundary_condition

  if boundary.type == "Neumann":  # TODO add support for more choice in Neumann
    mypb.setNeumannBoundaryCondition("GAUCHE")
    mypb.setNeumannBoundaryCondition("DROITE")
    mypb.setNeumannBoundaryCondition("HAUT")
    mypb.setNeumannBoundaryCondition("BAS")
  else:
    raise Exception("TODO boundary type unsuported %s" % boundary.type)

  numerical = diffusion.numerical_parameters

  mypb.setTimeScheme(getMethod(numerical))
  mypb.setLinearSolver(getLinear(numerical), getPreconditioner(numerical), getMaxIter(numerical))

  computation = diffusion.computation_parameters
  mypb.setCFL(computation.cfl)
  mypb.setPrecision(computation.precision)
  mypb.setMaxNbOfTimeStep(computation.max_time_step)
  mypb.setTimeMax(computation.max_time)
  mypb.setFreqSave(computation.save_frequency)
  mypb.setFileName(computation.file_name)
  #  TODO mypb.setResultDirectory(result_directory)
  mypb.setSaveFileFormat(getFormat(computation))

  # evolution
  mypb.initialize()
  logger.info("Run the simulation")
  ok = mypb.run()
  mypb.terminate()

  if ok:
    logger.info("Simulation completed status: Success")
  else:
    logger.info("Simulation completed status: Failed")


  return ok
