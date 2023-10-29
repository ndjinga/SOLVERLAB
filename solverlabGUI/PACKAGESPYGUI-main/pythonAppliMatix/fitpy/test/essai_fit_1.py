#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
see https://pypi.org/project/python-fit/
python 2-3 compliant
"""

import fitpy.fit as FIT
import numpy as np
import matplotlib.pyplot as plt

verbose = False

def essai_logXSurX(x, n):
    crossing = 2.
    xx = x + 1.
    return np.log(xx)/((xx/(crossing + 1.))**n)

def essai_plot():
  x = np.linspace( 0, 5, 1000)
  for n in [1, 2, 3, 4]:
    y = essai_logXSurX(x,n)
    plt.plot(x, y, 'b-', label='n=%s' % n)
  plt.legend()
  plt.show()
  return

if verbose: essai_plot()


def essai_1(isplot=True):
  # Create some data to fit
  x = np.arange(-10, 10, .2)
  hh = 1.
  cc = 2.
  ww = 3.
  nn = 2./10.
  # A gaussian of height hh, centered at cc, width ww, centered at cc, with noise nn.
  # y = hh*np.exp(-(x-cc)**ww/8) + (np.random.rand(100) - 0.5)*nn
  y = hh*np.exp(-0.5 * ((x-cc)/ww)**2) + (np.random.rand(100) - 0.5)*nn

  # No need to provide first guess at parameters for fit.gaus
  res = FIT.fit(FIT.gaus, x, y)
  (xf, yf), params, err, chi = res.getValues()

  if isplot:
    import matplotlib.pyplot as plt

    print("\ngaussian of height %s, centered at %s, width %s, with noise %s" % (hh, cc, ww, nn))
    print("height:    %5.2f +/- %5.3f" % (params[0], err[0]))
    print("centered:  %5.2f +/- %5.3f" % (params[1], err[1]))
    print("width:     %5.2f +/- %5.3f" % (params[2], err[2]))

    plt.plot(x, y, 'bo', label='Data')
    plt.plot(xf, yf, 'r-', label='Fit')
    plt.legend()
    plt.show()

  return res

def essai_2(isplot=True):
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
  y = np.array(y) + (np.random.rand(len(y)) - 0.5)*.05 # with noise

  # Gaussian simple
  def example_function_gaus(params, x, firstEstimate=None):
    """
    An un-normalized Gaussian curve.
    params: N, mu, sigma
    """
    if firstEstimate is None:
      N, mu, sigma = params
      return N * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    else:
      xdata = firstEstimate[0]
      ydata = firstEstimate[1]
      # The current behaviour of 'Series.argmax' is deprecated, use 'idxmax' instead.
      # This situation has occurred in the case of a downstream library like 'pandas'
      idx = np.argmax(np.array(ydata)) # downstream library like 'pandas'
      mu = xdata[idx]
      # print("idx %s mu %s\n%s" % (idx, mu, PP.pformat(dir(ydata))))
      # print("example_function_gaus idx %s mu %s" % (idx, mu))
      sigma = xdata.std()
      N = max(ydata)
      return np.array([N, mu, sigma])

  # Gaussian with quadratic background.
  def example_function_gaussian_quadratic(params, x, firstEstimate=None):
    """gaussian with quadratic background."""
    if firstEstimate is None:
      hh, cc, ww, a, b, c = params
      return hh * np.exp(-0.5 * ((x - cc) / ww) ** 2) + a * x ** 2 + b * x + c
    else:  # firstEstimate
      xdata = firstEstimate[0]
      ydata = firstEstimate[1]
      # The current behaviour of 'Series.argmax' is deprecated, use 'idxmax' instead.
      idx = np.argmax(np.array(ydata)) # downstream library like 'pandas'
      mu = xdata[idx]
      # print("idx %s mu %s\n%s" % (idx, mu, PP.pformat(dir(ydata))))
      # print("example_function_gaussian_quadratic idx %s mu %s" % (idx, mu))
      sigma = xdata.std()
      N = max(ydata)
      return np.array([N, mu, sigma, 0., 0., 0.])


  # It will still try to guess parameters, but they are dumb!
  resq = FIT.fit(example_function_gaussian_quadratic, x, y)
  # It will still try to guess parameters, but they are dumb!
  ress = FIT.fit(example_function_gaus, x, y)

  if isplot:
    import matplotlib.pyplot as plt
    (xfq, yfq), params, err, chi = resq.getValues()
    print("\ngaussian with quadratic background")
    print("height:    %5.2f +/- %5.3f" % (params[0], err[0]))
    print("centered:  %5.2f +/- %5.3f" % (params[1], err[1]))
    print("width:     %5.2f +/- %5.3f" % (params[2], err[2]))
    print("a:         %5.2f +/- %5.3f" % (params[3], err[0]))
    print("b:         %5.2f +/- %5.3f" % (params[4], err[1]))
    print("c:         %5.2f +/- %5.3f" % (params[5], err[2]))

  (xfs, yfs), params, err, chi = ress.getValues()

  if isplot:
    (xfs, yfs), params, err, chi = ress
    print("\ngaussian simple")
    print("height:    %5.2f +/- %5.3f" % (params[0], err[0]))
    print("centered:  %5.2f +/- %5.3f" % (params[1], err[1]))
    print("width:     %5.2f +/- %5.3f" % (params[2], err[2]))

    plt.plot(x, y, 'bo', label='Data')
    plt.plot(xfq, yfq, 'r-', label='Fit_quadra')
    plt.plot(xfs, yfs, 'b-', label='Fit_simple')
    plt.legend()
    plt.show()

  return resq

def essai_boltzmann():
  from matplotlib import pyplot as plt
  import numpy as np
  f = np.zeros(40)
  v = np.arange(0,4,0.1)
  for i in np.arange(0, 40):
      f[i] = v[i]**2*(np.exp(-v[i]**2))
  plt.plot(v,f)
  plt.show()

def essai_gamma():
  # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.gamma.html
  from scipy.stats import gamma
  from matplotlib import pyplot as plt
  import numpy as np
  f = np.zeros(40)
  v = np.arange(0,4,0.1)
  rv = gamma(3., loc = 0., scale = 2.)
  # produces a frozen form of gamma with shape a = 3., loc =0. and lambda = 1./scale = 1./2..
  print(rv)
  #
  a = 1.99323054838
  mean, var, skew, kurt = gamma.stats(a, moments='mvsk')
  # Display the probability density function (pdf):
  x = np.linspace(gamma.ppf(0.01, a), gamma.ppf(0.99, a), 30)
  print(x)
  plt.plot(x, gamma.pdf(x, a), 'o-', lw=5, alpha=0.6, label='gamma pdf')
  plt.show()


if __name__ == "__main__":
  np.random.seed(0)
  essai_1()
  essai_2()
  essai_plot()
