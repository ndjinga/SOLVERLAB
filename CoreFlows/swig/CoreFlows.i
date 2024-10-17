%module CoreFlows

%include std_string.i
%include std_vector.i
%include std_map.i

%include "Mesh.hxx" /* To include the enum EntityType */

#ifdef MPI4PY_ROOT_DIR
  %include mpi4py.i
  %mpi4py_typemap(Comm, MPI_Comm);
#endif

%apply bool& INOUT {bool &stop}

namespace std {
 %template(VectorDouble) vector<double>;
 %template(VectorVectorDouble) vector< vector<double> >;
 %template(VectorInt) vector<int>;
 %template(MapIntDouble) map<int,double>;
};

%{

#include "DriftModel.hxx"
#include "FiveEqsTwoFluid.hxx"
#include "IsothermalTwoFluid.hxx"
#include "ProblemFluid.hxx"
#include "ProblemCoreFlows.hxx"
#include "TransportEquation.hxx"
#include "DiffusionEquation.hxx"
#include "StationaryDiffusionEquation.hxx"
#include "IsothermalSinglePhase.hxx"
#include "SinglePhase.hxx"
#include "Fluide.h"
#include "StiffenedGas.hxx"
#include "WaveStaggered.hxx"
#include "EulerBarotropicStaggered.hxx"


%}

%include "ProblemCoreFlows.hxx"
%include "ProblemFluid.hxx"
%include "DriftModel.hxx"
%include "FiveEqsTwoFluid.hxx"
%include "IsothermalTwoFluid.hxx"
%include "TransportEquation.hxx"
%include "DiffusionEquation.hxx"
%include "StationaryDiffusionEquation.hxx"
%include "IsothermalSinglePhase.hxx"
%include "SinglePhase.hxx"
%include "Fluide.h"
%include "StiffenedGas.hxx"
%include "WaveStaggered.hxx"
%include "EulerBarotropicStaggered.hxx"

