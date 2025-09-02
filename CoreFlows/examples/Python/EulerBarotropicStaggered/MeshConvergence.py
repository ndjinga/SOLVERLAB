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

names = { "RTExplicit"  : pream.names["dSh_Quad_DG0_Godunov"],
          "RTImpl"   : pream.names["dSh_Quad_DG0_LaxFriedrich"],
          "RTHodgeLapla"    : pream.names["dQ_Quad_DG0_LaxFriedrich"],
          "dQ_Mach1em4_Quad_DG0_RoeBarotropic"   : pream.names["dQ_Quad_DG0_Godunov"]}
import math
from  matplotlib import pyplot as plt
import pandas as pd
import os 
import numpy as np
from addTriangleSlope import *

if __name__ == "__main__":
	pyplot.clf()
	# Wave Explicit GradDiv 
	NormL2cellsx = [0.37113782820166924, 0.19295793604330752,0.06602506462801658, 0.01926450110921030]
	NormL2cellsy = [0.12454330601332075, 0.07058931544256057,0.02862230869389324, 0.00909828393739803]
	
	# Wave Implicit  
	NormL2cellsx2 = [0.37113792871438178, 0.19295793813180534,0.06602509058398934, 0.01926450160227236]
	NormL2cellsy2 = [0.12454337314874667, 0.07058931742818185,0.02862234075815582, 0.00909828647447488]	
	
	# Wave Explicit HodgeLapla  
	NormL2cellsx3 = [3.23134478777776879, 3.08182233219911250,3.17571342589987982, 3.34127008684455706]
	NormL2cellsy3 = [1.27987038901922112, 1.43618120909834723,1.72899911734330303, 1.90073413571409300]	
	
	
	sizeMesh = [1.0/40, 1.0/80, 1.0/320, 1.0/1280]
	
	
	plt.figure()
	ylabel=r"$\Vert (u_y)_h^{\mathrm{Hodge-Lapla}} - (u_{y})^{\mathrm{ex}}\Vert_2$"
	nameFigure="_ConvergenceU.pdf"
	locs = "best"
	identity = "RTExplicit"
	identity2 = "RTImpl"
	identity3 = "RTHodgeLapla"
	
	pyplot.loglog(sizeMesh, NormL2cellsy3, c=names[identity3]["color"], ls=names[identity3]["ls"], marker=names[identity3]["marker"], markersize= names[identity3]["markersize"], 			markerfacecolor=names[identity3]["markerfacecolor"], markevery= 1, 	label=names[identity3]["flux"]+", "+names[identity3]["EF"])
	#pyplot.loglog(sizeMesh, NormL2cellsy, c=names[identity]["color"], ls=names[identity]["ls"], marker=names[identity]["marker"], markersize= names[identity]["markersize"], 			markerfacecolor=names[identity]["markerfacecolor"], markevery= 1, 	label=names[identity]["flux"]+", "+names[identity]["EF"])
	#addBottomTriangleSlope(sizeMesh,NormL2cellsx3,0.1,0.1,pream.fontsizeSlope)
	
	pyplot.rc('xtick', labelsize=pream.fontsizeAnnotate)
	pyplot.rc('ytick', labelsize=pream.fontsizeAnnotate)
	#pyplot.xticks([0.1,0.7])
	pyplot.xlabel("h",fontsize=pream.fontsizeLabel,labelpad=1)          
	pyplot.ylabel(ylabel,fontsize=pream.fontsizeLabel,labelpad=1)     
	pyplot.legend(loc=locs,fontsize=pream.fontsizeLeg)

	nameFigure="ErrorYHodgeLapla_WaveStaggered_Cylinder_"+"Quad"+nameFigure
	pp = PdfPages(nameFigure)
	pyplot.savefig(pp, format='pdf',bbox_inches="tight")
	pp.close()
	
