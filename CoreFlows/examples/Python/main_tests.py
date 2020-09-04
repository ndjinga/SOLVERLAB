#!/usr/bin/env python
# -*-coding:utf-8 -*

import DiffusionEquation_1DHeatedRod
import DriftModel_1DRiemannProblem
import DriftModel_1DPressureLoss
import DriftModel_1DBoilingChannel
import DriftModel_1DBoilingAssembly
import DriftModel_1DDepressurisation
import DriftModel_1DPorosityJump
import DriftModel_1DVidangeReservoir
import DriftModel_2DPressureLoss
import DriftModel_2DInclinedBoilingChannel
import DriftModel_2DBoilingChannelBarrier
import DriftModel_2DInclinedBoilingChannelBarrier
import DriftModel_2DVidangeReservoir
import DriftModel_2DVidangeReservoirUnstructured
import DriftModel_2BranchesBoilingChannels
import DriftModel_3DBoilingChannelBarrier
import FiveEqsTwoFluid_1DBoilingChannel
import FiveEqsTwoFluid_1DBoilingAssembly
import FiveEqsTwoFluid_1DVidangeReservoir
import FiveEqsTwoFluid_2DInclinedBoilingChannel
import FiveEqsTwoFluid_2DInclinedSedimentation
import FiveEqsTwoFluid_2DVidangeReservoir
import IsothermalTwoFluid_1DSedimentation
import IsothermalTwoFluid_1DVidangeReservoir
import IsothermalTwoFluid_2DVidangeReservoir
import SinglePhase_1DRiemannProblem
import SinglePhase_1DDepressurisation
import SinglePhase_1DWaterHammer
import SinglePhase_1DHeatedChannel
import SinglePhase_1DHeatedAssembly
import SinglePhase_2DHeatedChannelInclined
import SinglePhase_2DWallHeatedChannel_ChangeSect
import SinglePhase_2DLidDrivenCavity
import SinglePhase_2DLidDrivenCavity_unstructured
import SinglePhase_2DSphericalExplosion_unstructured
import SinglePhase_2DVidangeReservoir
import SinglePhase_3DHeatDrivenCavity
import SinglePhase_2DThermalDiffusion
import SinglePhase_2BranchesHeatedChannels
import TransportEquation_1DHeatedChannel

def main():
		
	if (DiffusionEquation_1DHeatedRod.DiffusionEquation_1DHeatedRod()):
		print( "Simulation python " + "DiffusionEquation_1DHeatedRod" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DiffusionEquation_1DHeatedRod" + "  failed ! " );
		return 0
	if (DriftModel_1DRiemannProblem.DriftModel_1DRiemannProblem()):
		print( "Simulation python " + "DriftModel_1DRiemannProlem" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DRiemannproblem" + "  failed ! " );
		return 0
	if (DriftModel_1DPressureLoss.DriftModel_1DPressureLoss()):
		print( "Simulation python " + "DriftModel_1DPressureLoss" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DPressureLoss" + "  failed ! " );
		return 0
	if (DriftModel_1DBoilingChannel.DriftModel_1DBoilingChannel()):
		print( "Simulation python " + "DriftModel_1DBoilingChannel" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DBoilingChannel" + "  failed ! " );
		return 0
	if (DriftModel_1DBoilingAssembly.DriftModel_1DBoilingAssembly()):
		print( "Simulation python " + "DriftModel_1DBoilingAssembly" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DBoilingAssembly" + "  failed ! " );
		return 0
	if (DriftModel_1DDepressurisation.DriftModel_1DDepressurisation()):
		print( "Simulation python " + "DriftModel_1DDepressurisation" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DDepressurisation" + "  failed ! " );
		return 0

	if (DriftModel_1DPorosityJump.DriftModel_1DPorosityJump()):
		print( "Simulation python " + "DriftModel_1DPorosityJump" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DPorosityJump" + "  failed ! " );
		return 0

	if (DriftModel_1DVidangeReservoir.DriftModel_1DVidangeReservoir()):
		print( "Simulation python " + "DriftModel_1DVidangeReservoir" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_1DVidangeReservoir" + "  failed ! " );
		return 0

	if (DriftModel_2DPressureLoss.DriftModel_2DPressureLoss()):
		print( "Simulation python " + "DriftModel_2DPressureLoss" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DPressureLoss" + "  failed ! " );
		return 0

	if (DriftModel_2DInclinedBoilingChannel.DriftModel_2DInclinedBoilingChannel()):
		print( "Simulation python " + "DriftModel_2DInclinedBoilingChannel" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DInclinedBoilingChannel" + "  failed ! " );
		return 0

	if (DriftModel_2DBoilingChannelBarrier.DriftModel_2DBoilingChannelBarrier()):
		print( "Simulation python " + "DriftModel_2DBoilingChannelBarrier" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DBoilingChannelBarrier" + "  failed ! " );
		return 0

	if (DriftModel_2DInclinedBoilingChannelBarrier.DriftModel_2DInclinedBoilingChannelBarrier()):
		print( "Simulation python " + "DriftModel_2DInclinedBoilingChannelBarrier" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DInclinedBoilingChannelBarrier" + "  failed ! " );
		return 0

	if (DriftModel_2DVidangeReservoir.DriftModel_2DVidangeReservoir()):
		print( "Simulation python " + "DriftModel_2DVidangeReservoir" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DVidangeReservoir" + "  failed ! " );
		return 0

	if (DriftModel_2DVidangeReservoirUnstructured.DriftModel_2DVidangeReservoirUnstructured()):
		print( "Simulation python " + "DriftModel_2DVidangeReservoirUnstructured" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2DVidangeReservoirUnstructured" + "  failed ! " );
		return 0

	if (DriftModel_3DBoilingChannelBarrier.DriftModel_3DBoilingChannelBarrier()):
		print( "Simulation python " + "DriftModel_3DBoilingChannelBarrier" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_3DBoilingChannelBarrier" + "  failed ! " );
		return 0

	if (DriftModel_2BranchesBoilingChannels.DriftModel_2BranchesBoilingChannels()):
		print( "Simulation python " + "DriftModel_2BranchesBoilingChannels" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "DriftModel_2BranchesBoilingChannels" + "  failed ! " );
		return 0

	if (FiveEqsTwoFluid_1DBoilingChannel.FiveEqsTwoFluid_1DBoilingChannel()):
		print( "Simulation python " + "FiveEqsTwoFluid_1DBoilingChannel" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "FiveEqsTwoFluid_1DBoilingChannel" + "  failed ! " );
		return 0

	if (FiveEqsTwoFluid_1DBoilingAssembly.FiveEqsTwoFluid_1DBoilingAssembly()):
		print( "Simulation python " + "FiveEqsTwoFluid_1DBoilingAssembly" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "FiveEqsTwoFluid_1DBoilingAssembly" + "  failed ! " );
		return 0

#	if (FiveEqsTwoFluid_1DVidangeReservoir.FiveEqsTwoFluid_1DVidangeReservoir()):
#		print( "Simulation python " + "FiveEqsTwoFluid_1DVidangeReservoir" + " is successful !" );
#		pass
#	else:
#		print( "Simulation python " + "FiveEqsTwoFluid_1DVidangeReservoir" + "  failed ! " );
#		return 0

	if (FiveEqsTwoFluid_2DInclinedBoilingChannel.FiveEqsTwoFluid_2DInclinedBoilingChannel()):
		print( "Simulation python " + "FiveEqsTwoFluid_2DInclinedBoilingChanne" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "FiveEqsTwoFluid_2DInclinedBoilingChannel" + "  failed ! " );
		return 0

	if (FiveEqsTwoFluid_2DInclinedSedimentation.FiveEqsTwoFluid_2DInclinedSedimentation()):
		print( "Simulation python " + "FiveEqsTwoFluid_2DInclinedSedimentation" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "FiveEqsTwoFluid_2DInclinedSedimentation" + "  failed ! " );
		return 0

#	if (FiveEqsTwoFluid_2DVidangeReservoir.FiveEqsTwoFluid_2DVidangeReservoir()):
#		print( "Simulation python " + "FiveEqsTwoFluid_2DVidangeReservoir" + " is successful !" );
#		pass
#	else:
#		print( "Simulation python " + "FiveEqsTwoFluid_2DVidangeReservoir" + "  failed ! " );
#		return 0

	if (IsothermalTwoFluid_1DSedimentation.IsothermalTwoFluid_1DSedimentation()):
		print( "Simulation python " + "IsothermalTwoFluid_1DSedimentation" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "IsothermalTwoFluid_1DSedimentation" + "  failed ! " );
		return 0

	if (IsothermalTwoFluid_1DVidangeReservoir.IsothermalTwoFluid_1DVidangeReservoir()):
		print( "Simulation python " + "IsothermalTwoFluid_1DVidangeReservoir" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "IsothermalTwoFluid_1DVidangeReservoir" + "  failed ! " );
		return 0
	if (IsothermalTwoFluid_2DVidangeReservoir.IsothermalTwoFluid_2DVidangeReservoir()):
		print( "Simulation python " + "IsothermalTwoFluid_2DVidangeReservoir" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "IsothermalTwoFluid_2DVidangeReservoir" + "  failed ! " );
		return 0

	if (SinglePhase_1DRiemannProblem.SinglePhase_1DRiemannProblem()):
		print( "Simulation python " + "SinglePhase_1DRiemannProblem" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_1DRiemannProblem" + "  failed ! " );
		return 0

	if (SinglePhase_1DDepressurisation.SinglePhase_1DDepressurisation()):
		print( "Simulation python " + "SinglePhase_1DDepressurisation" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_1DDepressurisation" + "  failed ! " );
		return 0

	if (SinglePhase_1DWaterHammer.SinglePhase_1DWaterHammer()):
		print( "Simulation python " + "SinglePhase_1DWaterHammer" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_1DWaterHammer" + "  failed ! " );
		return 0

	if (SinglePhase_1DHeatedChannel.SinglePhase_1DHeatedChannel()):
		print( "Simulation python " + "SinglePhase_1DHeatedChannel" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_1DheatedChannel" + "  failed ! " );
		return 0

	if (SinglePhase_1DHeatedAssembly.SinglePhase_1DHeatedAssembly()):
		print( "Simulation python " + "SinglePhase_1DHeatedAssembly" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_1DHeatedAssembly" + "  failed ! " );
		return 0

	if (SinglePhase_2DWallHeatedChannel_ChangeSect.SinglePhase_2DWallHeatedChannel_ChangeSect()):
		print( "Simulation python " + "2DWallHeatedChannel_ChangeSect" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "2DWallHeatedChannel_ChangeSect" + "  failed ! " );
		return 0

	if (SinglePhase_2DHeatedChannelInclined.SinglePhase_2DHeatedChannelInclined()):
		print( "Simulation python " + "SinglePhase_2DHeatedChannelInclined" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DHeatedChannelInclined" + "  failed ! " );
		return 0

	if (SinglePhase_2DLidDrivenCavity.SinglePhase_2DLidDrivenCavity()):
		print( "Simulation python " + "SinglePhase_2DLidDrivenCavity" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DLidDrivenCavity" + "  failed ! " );
		return 0

	if (SinglePhase_2DLidDrivenCavity_unstructured.SinglePhase_2DLidDrivenCavity_unstructured()):
		print( "Simulation python " + "SinglePhase_2DLidDrivenCavity_unstructured" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DLidDrivenCavity_unstructured" + "  failed ! " );
		return 0

	if (SinglePhase_2DSphericalExplosion_unstructured.SinglePhase_2DSphericalExplosion_unstructured()):
		print( "Simulation python " + "SinglePhase_2DSphericalExplosion_unstructured" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DSphericalExplosion_unstructured" + "  failed ! " );
		return 0

	if (SinglePhase_2DVidangeReservoir.SinglePhase_2DVidangeReservoir()):
		print( "Simulation python " + "SinglePhase_2DVidangeReservoir" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DVidangeReservoir" + "  failed ! " );
		return 0

	if (SinglePhase_3DHeatDrivenCavity.SinglePhase_3DHeatDrivenCavity()):
		print( "Simulation python " + "SinglePhase_3DHeatDrivenCavity" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_3DHeatDrivenCavity" + "  failed ! " );
		return 0

	if (SinglePhase_2DThermalDiffusion.SinglePhase_2DThermalDiffusion()):
		print( "Simulation python " + "SinglePhase_2DThermalDiffusion" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2DThermalDiffusion" + "  failed ! " );
		return 0

	if (SinglePhase_2BranchesHeatedChannels.SinglePhase_2BranchesHeatedChannels()):
		print( "Simulation python " + "SinglePhase_2BranchesHeatedChannels" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "SinglePhase_2BranchesHeatedChannels" + "  failed ! " );
		return 0

	if (TransportEquation_1DHeatedChannel.TransportEquation_1DHeatedChannel()):
		print( "Simulation python " + "TransportEquation_1DHeatedChannel" + " is successful !" );
		pass
	else:
		print( "Simulation python " + "TransportEquation_1DHeatedChannel" + "  failed ! " );
		return 0

	print("All python tests successful")
	return 1

if __name__ == """__main__""":
    main()
