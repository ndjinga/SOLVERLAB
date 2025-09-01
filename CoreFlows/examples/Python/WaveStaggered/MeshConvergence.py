#!/usr/bin/env python
# -*-coding:utf-8 -*

from matplotlib import pyplot
from matplotlib import markers
from pylab import genfromtxt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

import preambule as pream
from addTriangleSlope import *


pyplot.rcParams.update(pream.params)

names = { "dSh_Mach1em4_Quad_DG0_RoeBarotropic"  : pream.names["dSh_Quad_DG0_Godunov"],
          "dSh_Mach1em4_Quad_DG0_LaxFriedrich"   : pream.names["dSh_Quad_DG0_LaxFriedrich"],
          "dQ_Mach1em4_Quad_DG0_LaxFriedrich"    : pream.names["dQ_Quad_DG0_LaxFriedrich"],
          "dQ_Mach1em4_Quad_DG0_RoeBarotropic"   : pream.names["dQ_Quad_DG0_Godunov"]}
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
from addTriangleSlope import *

if __name__ == "__main__":
	pyplot.clf()
	NormL2cellsx = [0.37113782820166924, 0.19295793604330752,0.06602506462801658, 0.01926450110921030]
	NormL2cellsy = [0.12454330601332075, 0.07058931544256057,0.02862230869389324, 0.00909828393739803]
	sizeMesh = [1.0/40, 1.0/80, 1.0/320, 1.0/1280]
	
	addBottomTriangleSlope(sizeMesh,NormL2cellsx,0.1,1.0,pream.fontsizeSlope)
	plt.figure()
	ylabel=r"$\Vert (u_x)_h - (u_{x})_{\mathrm{ex}}\Vert_2$"
	nameFigure="_ConvergenceU.pdf"
	locs = "best"
	identity = "dSh_Mach1em4_Quad_DG0_RoeBarotropic"

	pyplot.loglog(sizeMesh, NormL2cellsx, c=names[identity]["color"], ls=names[identity]["ls"], marker=names[identity]["marker"], markersize= names[identity]["markersize"], 			markerfacecolor=names[identity]["markerfacecolor"], markevery= 1, 	label=names[identity]["flux"]+", "+names[identity]["EF"])


	pyplot.rc('xtick', labelsize=pream.fontsizeAnnotate)
	pyplot.rc('ytick', labelsize=pream.fontsizeAnnotate)
	pyplot.xticks([0.1,0.7])
	pyplot.xlabel("h",fontsize=pream.fontsizeLabel,labelpad=1)          
	pyplot.ylabel(ylabel,fontsize=pream.fontsizeLabel,labelpad=1)     
	pyplot.legend(loc=locs,fontsize=pream.fontsizeLeg)

	nameFigure="Error_WaveStaggered_Cylinder_"+"Quad"+nameFigure
	pp = PdfPages(nameFigure)
	pyplot.savefig(pp, format='pdf',bbox_inches="tight")
	pp.close()
	
