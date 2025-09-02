import numpy as N
from matplotlib import pyplot

def addTopTriangleSlope(xx,yy,diffy,diffannotate,fontsize):
    n=len(xx)
    #print("xx :",xx)
    #print("yy :",yy)
    xmin = 1.1 * xx[n-1]
    xmax = 0.9 * xx[n-2]
    deltax = xmax - xmin;
    slope = (N.log(yy[n-2])-N.log(yy[n-1]))/((N.log(xx[n-2])-N.log(xx[n-1])))
    fmax = yy[n-2] * (xmax/xx[n-2])**slope
    ymax = fmax * N.exp(diffy) 
    slope = round(slope,1)
    ymin = ymax*(xmin/xmax)**slope
    deltay = ymax - ymin;
    pyplot.arrow(xmin,ymin,deltax,deltay,linestyle="-",width=0)
    pyplot.arrow(xmin,ymin,0.,deltay,linestyle="-",width=0)
    pyplot.arrow(xmin,ymax,deltax,0.,linestyle="-",width=0)
    pyplot.annotate(r"slope="+str(slope),(xmin,ymax*N.exp(diffannotate)),fontsize=fontsize)


def addBottomTriangleSlope(xx,yy,diffy,diffannotate,fontsize):
    n=len(xx)
    #print("xx :",xx)
    #print("yy :",yy)
    xmin = 1.1 * xx[n-1]
    xmax = 0.9 * xx[n-2]
    deltax = xmax - xmin;
    slope = (N.log(yy[n-2])-N.log(yy[n-1]))/((N.log(xx[n-2])-N.log(xx[n-1])))
    fmax = yy[n-2] * (xmax/xx[n-2])**slope
    ymax = fmax * N.exp(-diffy) 
    slope = round(slope,1)
    ymin = ymax*(xmin/xmax)**slope
    deltay = ymax - ymin;
    pyplot.arrow(xmin,ymin,deltax,deltay,linestyle="-",width=0)
    pyplot.arrow(xmin,ymin,deltax,0.,linestyle="-",width=0)
    pyplot.arrow(xmax,ymin,0.,deltay,linestyle="-",width=0)
    pyplot.annotate(r"slope="+str(slope),(xmin,ymin*N.exp(-diffannotate)),fontsize=fontsize)
