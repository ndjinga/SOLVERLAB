#include "SinglePhase_1DRiemannProblem.cxx"
#include "SinglePhase_1DHeatedChannel.cxx"
#include "SinglePhase_1DDepressurisation.cxx"
#include "SinglePhase_2DWallHeatedChannel.cxx"
#include "SinglePhase_2DWallHeatedChannel_ChangeSect.cxx"
#include "SinglePhase_2DHeatedChannelInclined.cxx"
#include "SinglePhase_HeatedWire_2Branches.cxx"
#include "SinglePhase_2DLidDrivenCavity.cxx"
#include "SinglePhase_2DLidDrivenCavity_unstructured.cxx"
#include "SinglePhase_2DHeatDrivenCavity.cxx"
#include "SinglePhase_2DHeatDrivenCavity_unstructured.cxx"
#include "SinglePhase_2DSphericalExplosion_unstructured.cxx"
#include "SinglePhase_3DHeatDrivenCavity.cxx"
#include "DriftModel_1DBoilingChannel.cxx"
#include "DriftModel_1DBoilingAssembly.cxx"
#include "DriftModel_1DRiemannProblem.cxx"
#include "DriftModel_1DPressureLoss.cxx"
#include "DriftModel_1DDepressurisation.cxx"
#include "DriftModel_2DInclinedBoilingChannel.cxx"
#include "DriftModel_3DCanalCloison.cxx"
#include "IsothermalTwoFluid_1DSedimentation.cxx"
#include "IsothermalTwoFluid_1DRiemannProblem.cxx"
#include "IsothermalTwoFluid_1DDepressurisation.cxx"
#include "IsothermalTwoFluid_2DInclinedSedimentation.cxx"
#include "IsothermalTwoFluid_2DVidangeReservoir.cxx"
#include "FiveEqsTwoFluid_1DBoilingChannel.cxx"
#include "FiveEqsTwoFluid_1DRiemannProblem.cxx"
#include "FiveEqsTwoFluid_1DDepressurisation.cxx"
#include "FiveEqsTwoFluid_2DInclinedBoilingChannel.cxx"
#include "FiveEqsTwoFluid_2DInclinedSedimentation.cxx"
#include "DiffusionEquation_1DHeatedRod.cxx"
#include "TransportEquation_1DHeatedChannel.cxx"
#include "CoupledTransportDiffusionEquations_1DHeatedChannel.cxx"

using namespace std;


int main(int argc,char **argv)
{
	 if(!SinglePhase_1DRiemannProblem())
		throw CdmathException("test SinglePhase_1DRiemannProblem failed");
	else if(!SinglePhase_1DHeatedChannel())
		throw CdmathException("test SinglePhase_1DHeatedChannel() failed");
	else if(!SinglePhase_1DDepressurisation())
		throw CdmathException("test SinglePhase_1DDepressurisation() failed");
	else if (!SinglePhase_2DLidDrivenCavity())
		throw CdmathException("test SinglePhase_2DLidDrivenCavity failed");
	else if (!SinglePhase_2DLidDrivenCavity_unstructured())
		throw CdmathException("test SinglePhase_2DLidDrivenCavity_unstructured failed");
	else if (!SinglePhase_2DHeatDrivenCavity())
		throw CdmathException("test SinglePhase_2DHeatDrivenCavity failed");
	else if (!SinglePhase_2DHeatDrivenCavity_unstructured())
		throw CdmathException("test SinglePhase_2DHeatDrivenCavity_unstructured failed");
	else if(!SinglePhase_2DWallHeatedChannel())
		throw CdmathException("test SinglePhase_2DWallHeatedChannel() failed");
	else if(!SinglePhase_2DWallHeatedChannel_ChangeSect())
		throw CdmathException("test SinglePhase_2DWallHeatedChannel_ChangeSect() failed");
	else if(!SinglePhase_2DHeatedChannelInclined())
		throw CdmathException("test SinglePhase_2DHeatedChannelInclined() failed");
	else if(!SinglePhase_HeatedWire_2Branches())
		throw CdmathException("test SinglePhase_HeatedWire_2Branches() failed");
	else if (!SinglePhase_2DSphericalExplosion_unstructured())
		throw CdmathException("test SinglePhase_2DSphericalExplosion_unstructured failed");
	else if (!SinglePhase_3DHeatDrivenCavity())
		throw CdmathException("test SinglePhase_3DHeatDrivenCavity failed");
	else if(!DriftModel_1DRiemannProblem())
		throw CdmathException("test DriftModel_1DRiemannProblem failed ");
	else if(!DriftModel_1DPressureLoss())
		throw CdmathException("test DriftModel_1DPressureLoss failed ");
	else if(!DriftModel_1DBoilingChannel())
		throw CdmathException("test DriftModel_1DBoilingChannel failed ");
	else if(!DriftModel_1DBoilingAssembly())
		throw CdmathException("test DriftModel_1DBoilingAssembly failed ");
	else if(!DriftModel_1DDepressurisation())
		throw CdmathException("test DriftModel_1DDepressurisation failed ");
	else if(!DriftModel_2DInclinedBoilingChannel())
		throw CdmathException("test DriftModel_2DInclinedBoilingChannel failed ");
	else if(!DriftModel_3DCanalCloison())
		throw CdmathException("test DriftModel_3DCanalCloison failed ");
	else if(!IsothermalTwoFluid_1DRiemannProblem())
		throw CdmathException("test IsothermalTwoFluid_1DRiemannProblem failed");
	else if(!IsothermalTwoFluid_1DSedimentation())
		throw CdmathException("test IsothermalTwoFluid_1DSedimentation failed");
	else if(!IsothermalTwoFluid_1DDepressurisation())
		throw CdmathException("test IsothermalTwoFluid_1DDepressurisation failed");
	else if(!IsothermalTwoFluid_2DInclinedSedimentation())
		throw CdmathException("test IsothermalTwoFluid_2DInclinedSedimentation failed");
	else if(!IsothermalTwoFluid_2DVidangeReservoir())
		throw CdmathException("test IsothermalTwoFluid_2DVidangeReservoir failed");
	else if(!FiveEqsTwoFluid_1DRiemannProblem())
		throw CdmathException("test FiveEqsTwoFluid_1DRiemannProblem failed");
	else if(!FiveEqsTwoFluid_1DBoilingChannel())
		throw CdmathException("test FiveEqsTwoFluid_1DBoilingChannel failed");
	else if(!FiveEqsTwoFluid_1DDepressurisation())
		throw CdmathException("test FiveEqsTwoFluid_1DDepressurisation failed");
	else if(!FiveEqsTwoFluid_2DInclinedBoilingChannel())
		throw CdmathException("test FiveEqsTwoFluid_2DInclinedBoilingChannel failed ");
	else if(!FiveEqsTwoFluid_2DInclinedSedimentation())
		throw CdmathException("test FiveEqsTwoFluid_2DInclinedSedimentation failed ");
	else if(!TransportEquation_1DHeatedChannel())
		throw CdmathException("test TransportEquation_1DHeatedChannel() failed");
	else if(!DiffusionEquation_1DHeatedRod())
		throw CdmathException("test DiffusionEquation_1DHeatedRod() failed");
	else if(!CoupledTransportDiffusionEquations_1DHeatedChannel())
		throw CdmathException("test CoupledTransportDiffusionEquations_1DHeatedChannel() failed");//choose correct physical parameters to obtain good physical results
	else
		cout<<"All C tests successful"<<endl;

	return 1;
}
