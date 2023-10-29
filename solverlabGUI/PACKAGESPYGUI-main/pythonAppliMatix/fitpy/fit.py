#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
fit.py

modified python2-python3 compliant
from python-fit-1.0.0
     https://pypi.org/project/python-fit

A fitting package for python. Designed for ease of use for fast results.

List of built in functions that will try to return intelligent default
parameters for you when you use them:
  gaus
  expo
  double_exp
  line
  crystal_ball
  crystal_ball_norm

Built in functions that do not provide intelligent default values:
  fofoi - First over first inverse: a/(x+b) + c*x+d


================== Example Usage ===================================

from pylab import *
import fit
ion()


x = (4.2105303, 5.2631601, 6.2405997, 7.5187997, 8.7218, 
     9.7744402, 10.676691, 11.65414, 12.63158, 13.83459, 
     14.887219, 16.015039, 17.06767, 18.270679, 19.24812, 
     20.300751, 21.50376, 23.157888, 25.789471, 28.345871, 
     30.601501, 33.458643, 39.022559, 46.015039, 48.270679)

y = (0.18942001, 0.2099, 0.23891001, 0.27816002, 0.31911, 
     0.35836001, 0.39932001, 0.43686003, 0.46416002, 0.49829001, 
     0.51536004, 0.52556, 0.51876995, 0.5, 0.47271, 
     0.44026, 0.39249001, 0.33106002, 0.24060, 0.17746, 
     0.13311001, 0.11262, 0.095566, 0.095566, 0.095566)


x = np.array(x)
y = np.array(y)

# Guassian with linear background.
def example_function(params, x):
    N,mu,sigma,a,b,c = params
    return N*np.exp(-0.5 * ((x-mu)/sigma)**2 ) + a*x**2 + b*x + c

#f,p,e,chi = fit.fit(fit.expo, x,y)
(xf,yf),p,e,chi = fit.fit(example_function, x,y)
plot(x,y, 'b.')
plot(xf,yf)

=====================================================================

see for more:
https://sites.google.com/a/ucsc.edu/krumholz/teaching-and-courses/ast119_w15/class-10
see for next:
https://fr.wikipedia.org/wiki/Test_de_Kolmogorov-Smirnov
    -> scipy.stats.kstest with Python

"""

from scipy.odr import odrpack as odr
from scipy.odr import models

import numpy as np # exp, zeros, linspace, array, diff, average, logical_and, argmax
import pprint as PP

verbose = False

####################################################
class ReturnFit(object):
  """
  useful class to store fit information
  and method to print resume
  and method to display matplotlib
  """
  def getValues(self):
    """
    returns 'classic' contents fit values as ([xfit,yfit], params, err, chi)
    """
    params = self.coeff
    err = self.error
    chi = self.sum_square
    return self._fit, params, err, chi

  def __repr__(self):
    """avoiding print arrays in self._fit and self._dataToFit"""
    aDict = dict([ (k, v) for k, v in self.__dict__.items() if "_" != k[0] ])
    aDict["__doc__"] = self._func.__doc__
    res = PP.pformat(aDict)
    return "ReturnFit(\n %s\n)" % res[1:-1]

  def getX(self):
    """origin X of Y=f(X)"""
    return self._dataToFit[0]

  def getY(self):
    """origin Y of Y=f(X)"""
    return self._dataToFit[1]

  def getXf(self):
    """fitted X of fitted Y=f(X)"""
    return self._fit[0]

  def getYf(self):
    """fitted Y of fitted Y=f(X)"""
    return self._fit[1]


####################################################
def gaus(params, x, firstEstimate=None):
  """
  An un-normalized Gaussian curve.
  params: N, mu, sigma
  """
  if firstEstimate is None:
    N, mu, sigma = params
    return N*np.exp(-0.5 * ((x-mu)/sigma)**2 )
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    # The current behaviour of 'Series.argmax' is deprecated, use 'idxmax' instead.
    idx = np.argmax(np.array(ydata)) # downstream library like 'pandas'
    mu = xdata[idx]
    # print("idx %s mu %s\n%s" % (idx, mu, PP.pformat(dir(ydata))))
    # print("gaus idx %s mu %s" % (idx, mu))
    sigma = xdata.std()
    N = max(ydata)
    return np.array([N, mu, sigma])


def expo(params, x):
  """
  An exponential curve. *NB* The constant is in the exponent! 
  params: const, slope
  """
  if firstEstimate is None:
    const, slope = params
    return np.exp(const + slope*x)
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    const = np.average(xdata)
    slope = np.average(np.diff(ydata)/xdata[1:])
    return np.array([const, slope])


'''
def double_exp(params, x):
  """
  Two exponential constants in one.
  params: const1, slope1, mu1, const2, slope2, mu2
  """
  if firstEstimate is None:
    exp1 = params[0]*np.exp(params[1]*x + params[2])
    exp2 = params[3]*np.exp(params[4]*x + params[5])
    return exp1*exp2
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    const1, slope1, mu1, const2, slope2, mu2
    np.array([?, ?, ?, ?, ?, ?])
'''

def line(params, x):
  """
  Just another name for pol1.
  """
  if firstEstimate is None:
    intercept, slope = params
    return slope*x + intercept
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    return np.array([0, 0])


def crystal_ball(params, x):
  """
  A Gaussian curve on one side and a power-law on the other side. Used in
  physics to model lossy processes.
  See http://en.wikipedia.org/wiki/Crystal_Ball_function
  Note that the definition used here differs slightly. At the time of this
  writing, the wiki article has some discrepancies in definitions/plots. This
  definition makes it easier to fit the function by using complex numbers
  and by negating any negative values for a and n.

  This version of the crystal ball is normalized by an additional parameter.
  params: N, a, n, xb, sig
  """
  if firstEstimate is None:
    x = x+0j # Prevent warnings...
    N, a, n, xb, sig = params
    if a < 0:
      a = -a
    if n < 0:
      n = -n
    aa = abs(a)
    A = (n/aa)**n * np.exp(- aa**2 / 2)
    B = n/aa - aa
    total = 0.*x
    total += ((x-xb)/sig  > -a) * N * np.exp(- (x-xb)**2/(2.*sig**2))
    total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)
    try:
      return total.real
    except:
      return total
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    N = ydata.max()
    idx = np.argmax(np.array(ydata)) # downstream library like 'pandas'
    xb = xdata[idx]
    sigma = sum(ydata > ydata.mean()*1.8)*1. / len(xdata) * (xdata.max() - xdata.min())
    n = 2.
    a = .5
    return np.array([N, a, n, xb, sigma])


def crystal_ball_norm(params, x):
  """
  A Gaussian curve on one side and a power-law on the other side. Used in
  physics to model lossy processes.
  See http://en.wikipedia.org/wiki/Crystal_Ball_function
  Note that the definition used here differs slightly. At the time of this
  writing, the wiki article has some discrepancies in definitions/plots. This
  definition makes it easier to fit the function by using complex numbers
  and by negating any negative values for a and n.

  This version of the crystal ball is normalized by an internal normalization
  process.
  params: a, n, xb, sig
  """
  if firstEstimate is None:
    x = x + 0j # complex type, Prevent warnings...
    a, n, xb, sig = params
    if a < 0:
      a = -a
    if n < 0:
      n = -n
    aa = abs(a)
    A = (n/aa)**n * np.exp(- aa**2 / 2)
    B = n/aa - aa
    C = n/aa / (n-1.) * np.exp(-aa**2/2.)
    D = sqrt(pi/2.) * (1. + erf(aa/sqrt(2.)))
    N = 1. / (sig * (C+D))
    total = 0.*x
    total += ((x-xb)/sig  > -a) * N * np.exp(- (x-xb)**2/(2.*sig**2))
    total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)
    try:
      return total.real
    except:
      return total
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    idx = np.argmax(np.array(ydata)) # downstream library like 'pandas'
    xb = xdata[idx]
    sigma = sum(ydata > ydata.mean()*1.8)*1. / len(xdata) * (xdata.max() - xdata.min())
    n = 2.
    a = .5
    return np.array([a, n, xb, sigma])


def pow_law(params, x):
  """
  Power law curve.
  params: base,exponent
  """
  if firstEstimate is None:
    base,exponent = params
    return base*(x**exponent)
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]
    base = 2.6
    exponent = -.15
    return np.array([base, exponent])


def fofoi(params, x):
  """ fofoi - First over first inverse: a/(x+b) + c*x+d
  This function is useful when exponential curves and power-law curves don't
  quite capture a dual slope function.
  params: a, b, c, d
  """
  if firstEstimate is None:
    a, b, c, d = params
    return a/(x+b) + c*x+d
  else:
    xdata = firstEstimate[0]
    ydata = firstEstimate[1]


def fit(func, x, y, default_pars=None, data_range=None, we=None, verbose=False, itmax=300):
  """
  The meat of the fitting package. See docs of fit.py for more details.
  Functions available are gaus and expo and more.

  Error implementation provided via an example by: Tiago, 20071114

  Performs a least squares fit to the data, with errors!
  Uses scipy odrpack, but for least squares.
  
  IN: (func, x, y, verbose, itmax)
     func         - A function that accepts input in the form:
                    func(params, x)
     x,y (arrays) - data to fit
     default_pars - Optional default parameters to start minimization at.
     data_range   - Fit a subrange of (x,y). Provide a tuple of the form
                    (x_min, x_max).
     we           - Weighting for data points as delivered to ODR. You
                    probably want it the same length as x.
     verbose      - can be 0,1,2 for different levels of output
                    (False or True are the same as 0 or 1)
     itmax (int)  - optional maximum number of iterations
     
  OUT: (fit, params, err)
     fit    -  xfit,yfit arrays that can immediately plotted to see the
              results of your fit. xf and yf are defined as: 
                  xfit = np.linspace( min(x), max(x), len(x)*10)
                  yfit = func(xfit)
     params - the coefficients of your fit in the order the function takes.
     err    - standard error (1-sigma) on the coefficients
  """

  # If this is a histogram output, correct it for the user.
  if len(x) == len(y) + 1:
    x = (x[1:] + x[:-1])/2.
  # Take a slice of data.
  if data_range:
    y = y[np.logical_and(x > data_range[0], x < data_range[1])]
    x = x[np.logical_and(x > data_range[0], x < data_range[1])]

  # http://www.scipy.org/doc/api_docs/SciPy.odr.odrpack.html
  # see models.py and use ready made models!!!!
  if default_pars != None:
    beta0 = np.array(default_pars)
  else:
    # beta0 = get_default_params(x, y, func)
    beta0 = func(None, None, firstEstimate=[x, y])
  model_func  = models.Model(func)
  
  mydata = odr.Data(x, y, we)
  myodr  = odr.ODR(mydata, model_func, maxit=itmax, beta0=beta0)

  # Set type of fit to least-squares:
  myodr.set_job(fit_type=2)
  if verbose == 2: myodr.set_iprint(final=2)
        
  fit = myodr.run() # a type scipy.odr.odrpack.Output

  if False: #verbose: # print details debug
    # print("fit is %s\n%s" % (type(fit), PP.pformat(dir(fit))))
    # print("fit stop reason is %s" % PP.pformat(fit.stopreason))
    print("fit is: %s" % type(fit)) ; fit.pprint()

  if fit.stopreason[0] == 'Iteration limit reached':
    print('WARNING : poly_lsq Iteration limit reached, result not reliable.')

  # Results and errors
  coeff = fit.beta
  err   = fit.sd_beta
  chi   = fit.sum_square

  # The resulting fit.
  xfit = np.linspace( min(x), max(x), len(x)*10)
  yfit = func(fit.beta, xfit)

  res = ReturnFit() # casting for fit class
  res._func = func
  res._function_name = getFunctionName(func)
  res.coeff = coeff
  res.error = err
  res.sum_square = chi
  res._dataToFit = np.array([x,y])
  res._fit = np.array([xfit,yfit])
  if verbose:
    print("fit is:\n%s" % res)
  return res

def getFunctionName(func):
  try:     # python2
    return func.func_name
  except:  # python3
    return func.__name__


###########################################################
# Some more function definitions.
###########################################################

# First, create a polynomial decorator.

def polN_dec(func):
  """
  Polynomial decorator
  """
  degree = int(getFunctionName(func).replace("pol",""))
  def polN(*args):
    params, x = args
    if len(params) != degree+1:
      raise ValueError
    return sum( [param*x**i for i,param in enumerate(params)] )
  polN.__doc__ = """\
  A polynomial function of degree {}. Parameters are provided in
  increasing order of degree, i.e. [0] + [1]*x + [2]*x^2 + ...
  Input: (params, x)
    params - An array of length {}
    x      - The data to evaluate given the previous parameters.
  """
  polN.__doc__ = polN.__doc__.format(degree, degree+1)
  polN.__name__ = "pol{}".format(degree)
  return polN


# Create polynomials from 0 to 20.
@polN_dec
def pol0(params,x): return

@polN_dec
def pol1(params,x): return

@polN_dec
def pol2(params,x): return

@polN_dec
def pol3(params,x): return

@polN_dec
def pol4(params,x): return

@polN_dec
def pol5(params,x): return

@polN_dec
def pol6(params,x): return

@polN_dec
def pol7(params,x): return

@polN_dec
def pol8(params,x): return

@polN_dec
def pol9(params,x): return

@polN_dec
def pol10(params,x): return

@polN_dec
def pol11(params,x): return

@polN_dec
def pol12(params,x): return

@polN_dec
def pol13(params,x): return

@polN_dec
def pol14(params,x): return

@polN_dec
def pol15(params,x): return

@polN_dec
def pol16(params,x): return

@polN_dec
def pol17(params,x): return

@polN_dec
def pol18(params,x): return

@polN_dec
def pol19(params,x): return

@polN_dec
def pol20(params,x): return

