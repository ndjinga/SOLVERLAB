#include "utilitaire_algebre.h"
#include <iostream>
using namespace std;

roots_polynoms::roots_polynoms(): smalno(1.2e-38), infin(3.4e38), 
				  eta(2.22e-16), base(2) 
{
}

int roots_polynoms::rpoly(double *op, int degree, double *zeror, double *zeroi) 
{
  double t,aa,bb,cc,*temp,factor,rot;
  double *pt;
  double lo,max,min,xx,yy,cosr,sinr,xxx,x,sc,bnd;
  double xm,ff,df,dx;
  int cnt,nz,i,j,jj,l,nm1,zerok;

  /*  The following statements set machine constants. */
  are = eta;
  mre = eta;
  lo = smalno/eta;
  /*  Initialization of constants for shift rotation. */        
  xx = sqrt(0.5);
  yy = -xx;
  rot = 94.0;
  rot *= 0.017453293;
  cosr = cos(rot);
  sinr = sin(rot);
  n = degree;
  /*  Algorithm fails of the leading coefficient is zero. */
  if (fabs(op[0]) <=  smalno) return -1;
  /*  Remove the zeros at the origin, if any. */
  while (op[n] == 0.0) {
    j = degree - n;
    zeror[j] = 0.0;
    zeroi[j] = 0.0;
    n--;
  }
  if (n < 1) return -1;

  /*
   *  Allocate memory here
   */
  temp = new double [degree+1];
  pt = new double [degree+1];
  p = new double [degree+1];
  qp = new double [degree+1];
  k = new double [degree+1];
  qk = new double [degree+1];
  svk = new double [degree+1];
  /*  Make a copy of the coefficients. */
  for (i=0;i<=n;i++)
    p[i] = op[i];
  /*  Start the algorithm for one zero. */
 _40:        
  if (n == 1) {
    zeror[degree-1] = -p[1]/p[0];
    zeroi[degree-1] = 0.0;
    n -= 1;
    goto _99;
  }
  /*  Calculate the final zero or pair of zeros. */
  if (n == 2) {
    quad(p[0],p[1],p[2],&zeror[degree-2],&zeroi[degree-2],
	 &zeror[degree-1],&zeroi[degree-1]);
    n -= 2;
    goto _99;
  }
  /*  Find largest and smallest moduli of coefficients. */
  max = 0.0;
  min = infin;
  for (i=0;i<=n;i++) {
    x = fabs(p[i]);
    if (x > max) max = x;
    if (x != 0.0 && x < min){ min = x;
    //cerr << " x= " << x << " smalno " << smalno << endl;
    }
  }
  /*  Scale if there are large or very small coefficients.
   *  Computes a scale factor to multiply the coefficients of the
   *  polynomial. The scaling si done to avoid overflow and to
   *  avoid undetected underflow interfering with the convergence
   *  criterion. The factor is a power of the base.
   */
  sc = lo/min;
  if (sc > 1.0 && infin/sc < max) goto _110;
  if (sc <= 1.0) {
    if (max < 10.0) goto _110;
    if (sc == 0.0)
      sc = smalno;
  }
  l = (int)(log(sc)/log(base) + 0.5);
  factor = pow(base*1.0,l);
  if (factor != 1.0) {
    for (i=0;i<=n;i++) 
      p[i] = factor*p[i];     /* Scale polynomial. */
  }
 _110:
  /*  Compute lower bound on moduli of roots. */
  for (i=0;i<=n;i++) {
    pt[i] = (fabs(p[i]));
  }
  pt[n] = - pt[n];
  /*  Compute upper estimate of bound. */
  x = exp((log(-pt[n])-log(pt[0])) / (double)n);
  /*  If Newton step at the origin is better, use it. */        
  //if (pt[n-1]  !=  0.0) {
  if (fabs(pt[n-1] ) >  smalno*fabs(pt[n])) {
    xm = -pt[n]/pt[n-1];
    if (xm < x)  x = xm;
  } 
  /*  Chop the interval (0,x) until ff <= 0 */
  while (1) {
    xm = x*0.1;
    ff = pt[0];
    for (i=1;i<=n;i++) 
      ff = ff*xm + pt[i];
    if (ff <= 0.0) break;
    x = xm;
  }
  dx = x;
  /*  Do Newton interation until x converges to two 
   *  decimal places. 
   */
  while (fabs(dx/x) > 0.005) {
    ff = pt[0];
    df = ff;
    for (i=1;i<n;i++) { 
      ff = ff*x + pt[i];
      df = df*x + ff;
    }
    ff = ff*x + pt[n];
    dx = ff/df;
    x -= dx;
  }
  bnd = x;
  /*  Compute the derivative as the initial k polynomial
   *  and do 5 steps with no shift.
   */
  nm1 = n - 1;
  for (i=1;i<n;i++)
    k[i] = (double)(n-i)*p[i]/(double)n;
  k[0] = p[0];
  aa = p[n];
  bb = p[n-1];
  zerok = (fabs(k[n-1])< smalno);
  for(jj=0;jj<5;jj++) {
    cc = k[n-1];
    if (!zerok) {//cerr << "cc= " << cc << endl; 
      /*  Use a scaled form of recurrence if value of k at 0 is nonzero. */             
      t = -aa/cc;
      for (i=0;i<nm1;i++) {
	j = n-i-1;
	k[j] = t*k[j-1]+p[j];
      }
      k[0] = p[0];
      zerok = (fabs(k[n-1]) <= fabs(bb)*eta*10.0 || (fabs(k[n-1])< smalno));
    }
    else {
      /*  Use unscaled form of recurrence. */
      for (i=0;i<nm1;i++) {
	j = n-i-1;
	k[j] = k[j-1];
      }
      k[0] = 0.0;
      zerok = (fabs(k[n-1])< smalno);
    }
  }
  /*  Save k for restarts with new shifts. */
  for (i=0;i<n;i++) 
    temp[i] = k[i];
  /*  Loop to select the quadratic corresponding to each new shift. */
  for (cnt = 0;cnt < 20;cnt++) {
    /*  Quadratic corresponds to a double shift to a            
     *  non-real point and its complex conjugate. The point
     *  has modulus bnd and amplitude rotated by 94 degrees
     *  from the previous shift.
     */ 
    xxx = cosr*xx - sinr*yy;
    yy = sinr*xx + cosr*yy;
    xx = xxx;
    sr = bnd*xx;
    si = bnd*yy;
    u = -2.0 * sr;
    v = bnd;
    fxshfr(20*(cnt+1),&nz);
    if (nz != 0) {
      /*  The second stage jumps directly to one of the third
       *  stage iterations and returns here if successful.
       *  Deflate the polynomial, store the zero or zeros and
       *  return to the main algorithm.
       */
      j = degree - n;
      zeror[j] = szr;
      zeroi[j] = szi;
      n -= nz;
      for (i=0;i<=n;i++)
	p[i] = qp[i];
      if (nz != 1) {
	zeror[j+1] = lzr;
	zeroi[j+1] = lzi;
      }
      goto _40;
    }
    /*  If the iteration is unsuccessful another quadratic
     *  is chosen after restoring k.
     */
    for (i=0;i<n;i++) {
      k[i] = temp[i];
    }
  } 
  /*  Return with failure if no convergence with 20 shifts. */
 _99:
  delete [] svk;
  delete [] qk;
  delete [] k;
  delete [] qp;
  delete [] p;
  delete [] pt;
  delete [] temp;

//   cerr << " rpoly roots " << endl;
//   for (int ii=0; ii<degree;ii++) {
//     cerr << " (" << zeror[ii] << "," << zeroi[ii] << ")";
//   }
//   cerr << endl;
//   cerr << " rpoly roots point " << &zeror[0] << " " <<  &zeroi[0] << endl;
//   cerr << " rpoly roots  point 1 " << &zeror[1] << " " <<  &zeroi[1] << endl;


  return degree - n;
}

/*  Computes up to L2 fixed shift k-polynomials,
 *  testing for convergence in the linear or quadratic
 *  case. Initiates one of the variable shift
 *  iterations and returns with the number of zeros
 *  found.
 */
 void roots_polynoms::fxshfr(int l2,int *nz)
{
  double svu,svv,ui,vi,s;
  double betas,betav,oss,ovv,ss,vv,ts,tv;
  double ots,otv,tvv,tss;
  int type, i,j,iflag,vpass,spass,vtry,stry;

  *nz = 0;
  betav = 0.25;
  betas = 0.25;
  oss = sr;
  ovv = v;
  /*  Evaluate polynomial by synthetic division. */
  quadsd(n,&u,&v,p,qp,&a,&b);
  calcsc(&type);
  for (j=0;j<l2;j++) {
    /*  Calculate next k polynomial and estimate v. */
    nextk(&type);
    calcsc(&type);
    newest(type,&ui,&vi);
    vv = vi;
    /*  Estimate s. */
    ss = 0.0;
    if (k[n-1] != 0.0) ss = -p[n]/k[n-1];
    tv = 1.0;
    ts = 1.0;
    if (j == 0 || type == 3) goto _70;
    /*  Compute relative measures of convergence of s and v sequences. */
    if (vv != 0.0) tv = fabs((vv-ovv)/vv);
    if (ss != 0.0) ts = fabs((ss-oss)/ss);
    /*  If decreasing, multiply two most recent convergence measures. */
    tvv = 1.0;
    if (tv < otv) tvv = tv*otv;
    tss = 1.0;
    if (ts < ots) tss = ts*ots;
    /*  Compare with convergence criteria. */
    vpass = (tvv < betav);
    spass = (tss < betas);
    if (!(spass || vpass)) goto _70;
    /*  At least one sequence has passed the convergence test.
     *  Store variables before iterating.
     */
    svu = u;
    svv = v;
    for (i=0;i<n;i++) {
      svk[i] = k[i];
    }
    s = ss;
    /*  Choose iteration according to the fastest converging
     *  sequence.
     */
    vtry = 0;
    stry = 0;
    if (spass && (!vpass) || tss < tvv) goto _40;
  _20:        
    quadit(&ui,&vi,nz);
    if (*nz > 0) return;
    /*  Quadratic iteration has failed. Flag that it has
     *  been tried and decrease the convergence criterion.
     */
    vtry = 1;
    betav *= 0.25;
    /*  Try linear iteration if it has not been tried and
     *  the S sequence is converging.
     */
    if (stry || !spass) goto _50;
    for (i=0;i<n;i++) {
      k[i] = svk[i];
    }
  _40:
    realit(&s,nz,&iflag);
    if (*nz > 0) return;
    /*  Linear iteration has failed. Flag that it has been
     *  tried and decrease the convergence criterion.
     */
    stry = 1;
    betas *=0.25;
    if (iflag == 0) goto _50;
    /*  If linear iteration signals an almost double real
     *  zero attempt quadratic iteration.
     */
    ui = -(s+s);
    vi = s*s;
    goto _20;
    /*  Restore variables. */
  _50:
    u = svu;
    v = svv;
    for (i=0;i<n;i++) {
      k[i] = svk[i];
    }
    /*  Try quadratic iteration if it has not been tried
     *  and the V sequence is convergin.
     */
    if (vpass && !vtry) goto _20;
    /*  Recompute QP and scalar values to continue the
     *  second stage.
     */
    quadsd(n,&u,&v,p,qp,&a,&b);
    calcsc(&type);
  _70:
    ovv = vv;
    oss = ss;
    otv = tv;
    ots = ts;
  }
}
/*  Variable-shift k-polynomial iteration for a
 *  quadratic factor converges only if the zeros are
 *  equimodular or nearly so.
 *  uu, vv - coefficients of starting quadratic.
 *  nz - number of zeros found.
 */
 void roots_polynoms::quadit(double *uu,double *vv,int *nz)
{
  double ui,vi;
  double mp,omp,ee,relstp,t,zm;
  int type,i,j,tried;

  *nz = 0;
  tried = 0;
  u = *uu;
  v = *vv;
  j = 0;
  /*  Main loop. */
 _10:    
  quad(1.0,u,v,&szr,&szi,&lzr,&lzi);
  /*  Return if roots of the quadratic are real and not
   *  close to multiple or nearly equal and of opposite
   *  sign.
   */
  if (fabs(fabs(szr)-fabs(lzr)) > 0.01 * fabs(lzr)) return;
  /*  Evaluate polynomial by quadratic synthetic division. */
  quadsd(n,&u,&v,p,qp,&a,&b);
  mp = fabs(a-szr*b) + fabs(szi*b);
  /*  Compute a rigorous bound on the rounding error in
   *  evaluating p.
   */
  zm = sqrt(fabs(v));
  ee = 2.0*fabs(qp[0]);
  t = -szr*b;
  for (i=1;i<n;i++) {
    ee = ee*zm + fabs(qp[i]);
  }
  ee = ee*zm + fabs(a+t);
  ee *= (5.0 *mre + 4.0*are);
  ee = ee - (5.0*mre+2.0*are)*(fabs(a+t)+fabs(b)*zm)+2.0*are*fabs(t);
  /*  Iteration has converged sufficiently if the
   *  polynomial value is less than 20 times this bound.
   */
  if (mp <= 20.0*ee) {
    *nz = 2;
    return;
  }
  j++;
  /*  Stop iteration after 20 steps. */
  if (j > 20) return;
  if (j < 2) goto _50;
  if (relstp > 0.01 || mp < omp || tried) goto _50;
  /*  A cluster appears to be stalling the convergence.
   *  Five fixed shift steps are taken with a u,v close
   *  to the cluster.
   */
  if (relstp < eta) relstp = eta;
  relstp = sqrt(relstp);
  u = u - u*relstp;
  v = v + v*relstp;
  quadsd(n,&u,&v,p,qp,&a,&b);
  for (i=0;i<5;i++) {
    calcsc(&type);
    nextk(&type);
  }
  tried = 1;
  j = 0;
 _50:
  omp = mp;
  /*  Calculate next k polynomial and new u and v. */
  calcsc(&type);
  nextk(&type);
  calcsc(&type);
  newest(type,&ui,&vi);
  /*  If vi is zero the iteration is not converging. */
  if (vi == 0.0) return;
  relstp = fabs((vi-v)/vi);
  u = ui;
  v = vi;
  goto _10;
}
/*  Variable-shift H polynomial iteration for a real zero.
 *  sss - starting iterate
 *  nz  - number of zeros found
 *  iflag - flag to indicate a pair of zeros near real axis.
 */
 void roots_polynoms::realit(double *sss, int *nz, int *iflag)
{
  double pv,kv,t,s;
  double ms,mp,omp,ee;
  int i,j;

  *nz = 0;
  s = *sss;
  *iflag = 0;
  j = 0;
  /*  Main loop */
  while (1) {
    pv = p[0];
    /*  Evaluate p at s. */
    qp[0] = pv;
    for (i=1;i<=n;i++) {
      pv = pv*s + p[i];
      qp[i] = pv;
    }
    mp = fabs(pv);
    /*  Compute a rigorous bound on the error in evaluating p. */
    ms = fabs(s);
    ee = (mre/(are+mre))*fabs(qp[0]);
    for (i=1;i<=n;i++) {
      ee = ee*ms + fabs(qp[i]);
    }
    /*  Iteration has converged sufficiently if the polynomial
     *  value is less than 20 times this bound.
     */
    if (mp <= 20.0*((are+mre)*ee-mre*mp)) {
      *nz = 1;
      szr = s;
      szi = 0.0;
      return;
    }
    j++;
    /*  Stop iteration after 10 steps. */
    if (j > 10) return;
    if (j < 2) goto _50;
    if (fabs(t) > 0.001*fabs(s-t) || mp < omp) goto _50;
    /*  A cluster of zeros near the real axis has been
     *  encountered. Return with iflag set to initiate a
     *  quadratic iteration.
     */
    *iflag = 1;
    *sss = s;
    return;
    /*  Return if the polynomial value has increased significantly. */
  _50:
    omp = mp;
    /*  Compute t, the next polynomial, and the new iterate. */
    kv = k[0];
    qk[0] = kv;
    for (i=1;i<n;i++) {
      kv = kv*s + k[i];
      qk[i] = kv;
    }
    if (fabs(kv) <= fabs(k[n-1])*10.0*eta) {
      /*  Use unscaled form. */
      k[0] = 0.0;
      for (i=1;i<n;i++) {
	k[i] = qk[i-1];
      }
    }
    else {
      /*  Use the scaled form of the recurrence if the value
       *  of k at s is nonzero.
       */
      t = -pv/kv;
      k[0] = qp[0];
      for (i=1;i<n;i++) {
	k[i] = t*qk[i-1] + qp[i];
      }
    }
    kv = k[0];
    for (i=1;i<n;i++) {
      kv = kv*s + k[i];
    }
    t = 0.0;
    if (fabs(kv) > fabs(k[n-1]*10.0*eta)) t = -pv/kv;
    s += t;
  }
}

/*  This routine calculates scalar quantities used to
 *  compute the next k polynomial and new estimates of
 *  the quadratic coefficients.
 *  type - integer variable set here indicating how the
 *  calculations are normalized to avoid overflow.
 */
 void roots_polynoms::calcsc(int *type)
{
  /*  Synthetic division of k by the quadratic 1,u,v */    
  quadsd(n-1,&u,&v,k,qk,&c,&d);
  if (fabs(c) > fabs(k[n-1]*100.0*eta)) goto _10;
  if (fabs(d) > fabs(k[n-2]*100.0*eta)) goto _10;
  *type = 3;
  /*  Type=3 indicates the quadratic is almost a factor of k. */
  return;
 _10:
  if (fabs(d) < fabs(c)) {
    *type = 1;
    /*  Type=1 indicates that all formulas are divided by c. */   
    e = a/c;
    f = d/c;
    g = u*e;
    h = v*b;
    a3 = a*e + (h/c+g)*b;
    a1 = b - a*(d/c);
    a7 = a + g*d + h*f;
    return;
  }
  *type = 2;
  /*  Type=2 indicates that all formulas are divided by d. */
  e = a/d;
  f = c/d;
  g = u*b;
  h = v*b;
  a3 = (a+g)*e + h*(b/d);
  a1 = b*f-a;
  a7 = (f+u)*a + h;
}
/*  Computes the next k polynomials using scalars 
 *  computed in calcsc.
 */
 void roots_polynoms::nextk(int *type)
{
  double temp;
  int i;

  if (*type == 3) {
    /*  Use unscaled form of the recurrence if type is 3. */
    k[0] = 0.0;
    k[1] = 0.0;
    for (i=2;i<n;i++) {
      k[i] = qk[i-2];
    }
    return;
  }
  temp = a;
  if (*type == 1) temp = b;
  if (fabs(a1) <= fabs(temp)*eta*10.0) {
    /*  If a1 is nearly zero then use a special form of the
     *  recurrence.
     */
    k[0] = 0.0;
    k[1] = -a7*qp[0];
    for(i=2;i<n;i++) {
      k[i] = a3*qk[i-2] - a7*qp[i-1];
    }
    return;
  }
  /*  Use scaled form of the recurrence. */
  a7 /= a1;
  a3 /= a1;
  k[0] = qp[0];
  k[1] = qp[1] - a7*qp[0];
  for (i=2;i<n;i++) {
    k[i] = a3*qk[i-2] - a7*qp[i-1] + qp[i];
  }
}
/*  Compute new estimates of the quadratic coefficients
 *  using the scalars computed in calcsc.
 */
 void roots_polynoms::newest(int type,double *uu,double *vv)
{
  double a4,a5,b1,b2,c1,c2,c3,c4,temp;
  //  double smalno = 1.2e-38;
  //  double smalno = 0.0;

  /* Use formulas appropriate to setting of type. */
  if (type == 3) {
    /*  If type=3 the quadratic is zeroed. */
    *uu = 0.0;
    *vv = 0.0;
    return;
  }
  if (type == 2) {
    a4 = (a+g)*f + h;
    a5 = (f+u)*c + v*d;
  }
  else {
    a4 = a + u*b +h*f;
    a5 = c + (u+v*f)*d;
  }
  /*  Evaluate new quadratic coefficients. */
  b1 = -k[n-1]/p[n];
  b2 = -(k[n-2]+b1*p[n-1])/p[n];
  c1 = v*b2*a1;
  c2 = b1*a7;
  c3 = b1*b1*a3;
  c4 = c1 - c2 - c3;
  temp = a5 + b1*a4 - c4;
  if (fabs(temp) <= smalno) {
    *uu = 0.0;
    *vv = 0.0;
    return;
  }
  *uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7))/temp;
  *vv = v*(1.0+c4/temp);
  return;
}

/*  Divides p by the quadratic 1,u,v placing the quotient
 *  in q and the remainder in a,b.
 */
 void roots_polynoms::quadsd(int my_n,double *my_u,double *my_v,
					 double *my_p,double *my_q,
					 double *my_a,double *my_b)
{
  double my_c;
  int i;
  *my_b = my_p[0];
  my_q[0] = *my_b;
  *my_a = my_p[1] - (*my_b)*(*my_u);
  my_q[1] = *my_a;
  for (i=2;i<=my_n;i++) {
    my_c = my_p[i] - (*my_a)*(*my_u) - (*my_b)*(*my_v);
    my_q[i] = my_c;
    *my_b = *my_a;
    *my_a = my_c;
  }
}
/*  Calculate the zeros of the quadratic a*z^2 + b1*z + c.
 *  The quadratic formula, modified to avoid overflow, is used 
 *  to find the larger zero if the zeros are real and both
 *  are complex. The smaller real zero is found directly from 
 *  the product of the zeros c/a.
 */
 void roots_polynoms::quad(double my_a,double my_b1,double my_c,double *my_sr,double *my_si,
				       double *my_lr,double *my_li)
{
  double my_b,my_d,my_e;
  //  double smalno = 1.2e-38;

  if (fabs(my_a) <= smalno) {         /* less than two roots */
    if (fabs(my_b1) >= smalno)     
      *my_sr = -my_c/my_b1;
    else 
    *my_sr = 0.0;
    *my_lr = 0.0;
    *my_si = 0.0;
    *my_li = 0.0;
    return;
  }
  if (fabs(my_c) <= smalno) {         /* one real root, one zero root */
    *my_sr = 0.0;
    *my_lr = -my_b1/my_a;
    *my_si = 0.0;
    *my_li = 0.0;
    return;
  }
  /* Compute discriminant avoiding overflow. */
  my_b = my_b1/2.0;
  if (fabs(my_b) < fabs(my_c)) { 
    if (my_c < 0.0) 
      my_e = -my_a;
    else
      my_e = my_a;
    my_e = my_b*(my_b/fabs(my_c)) - my_e;
    my_d = sqrt(fabs(my_e))*sqrt(fabs(my_c));
  }
  else {
    my_e = 1.0 - (my_a/my_b)*(my_c/my_b);
    my_d = sqrt(fabs(my_e))*fabs(my_b);
  }
  if (my_e < 0.0) {      /* complex conjugate zeros */
    *my_sr = -my_b/my_a;
    *my_lr = *my_sr;
    *my_si = fabs(my_d/my_a);
    *my_li = -(*my_si);
  }
  else {
    if (my_b >= 0.0)   /* real zeros. */
      my_d = -my_d;
    *my_lr = (-my_b+my_d)/my_a;
    *my_sr = 0.0;
    if (*my_lr != 0.0) 
      *my_sr = (my_c/ *my_lr)/my_a;
    *my_si = 0.0;
    *my_li = 0.0;
  }
}

 //qques outils d'algebre lineaire

 /** \fn add
  * \brief Calcule la somme de 2vecteurs ( U=U+V)
  * @param *U,*V les 2 vecteurs
  * @param N taille des 2 vecteurs
  * @return la somme de U et V dans U*/
 void Polynoms::add
 (
 		double *U,		//vecteur auquel on ajoute
 		int N,			//taille du vecteur
 		const double *V	//ajout
 )
 {
 	for(int i=0;i<N;i++) U[i] += V[i];
 };

 /** \fn tensor
  * \brief Calcule le tenseur de 2 vecteurs
  * @param *a le 1er vecteur (*double)
  * @param N la taille du premier vecteur (int)
  * @param *b le 2ieme vecteur (*double)
  * @param M la taille du 2ieme vecteur (int)
  * @param *ab le resultat
  * @return */

 void Polynoms::tensor
 (	const double *a,	//vecteur gauche
 		int N,			//taille du vecteur gauche
 		const double *b,	//vecteur droit
 		int M,			//taille du vecteur droit
 		double *ab		//conteneur de sortie
 )
 {
 	for(int i=0;i<N;i++)
 		for(int j=0;j<M;j++)
 			ab[i*M+j] = a[i]*b[j];
 };

 /** \fn matrixProduct
  * \brief calcule le produit de 2 matrices .
  * @param *Matrix1 la 1ere matrice
  * @param R1,C1 la taille de la matrice1 Nbr de ligne et le nbr de colonne
  * @param Matrix2 la 2ieme matrice
  * @param R2,C2 la taille de la matrice2 Nbr de ligne et le nbr de colonne
  * @param out le resultat de Matrix1*Matrix2
  * @return */
 void Polynoms::matrixProduct
 (
 		const double *Matrix1,
 		const int &R1,
 		const int &C1,
 		const double *Matrix2,
 		const int &R2,
 		const int &C2,
 		double *out
 )
 {
 	for(int i=0; i<R1; i++)
 	{
 		for(int j=0; j<C2; j++)
 		{
 			out[i*C2+j] = 0;
 			for(int k=0; k<C1; k++)
 			{
 				out[i*C2+j] += Matrix1[i*C1+k]*Matrix2[k*C2+j];
 			}
 		}
 	}
 };
 /** \fn matrixProdVec
  * \brief Calcule le produit Matrice vecteur
  * @param *Matrix correspond à la matrice du produit
  * @param R1,C1 la taille de la matrice( Nbr de ligne et le nbr de colonne)
  * @param Vector correspond au vecteur du produit
  * @param out le resultat de Matrix*Vector
  * @return   */
 void Polynoms::matrixProdVec
 (	const double *Matrix,
 		const int &R1,
 		const int &C1,
 		const double *Vector,
 		double *out
 )
 {
 	for(int i=0; i<R1; i++)
 	{
 		out[i] = 0;
 		for(int j=0; j<C1; j++)
 		{
 			out[i] += Matrix[i*C1+j]*Vector[j];
 		}
 	}
 };

 /** \fn shift_diagonal
  * \brief Calcule la trace d'une matrice carrée
  * @param *Matrix correspond  à la matrice
  * @param N taille de la matrice
  * @param shift résultat
  * @return renvoie la trace de la matrice Matrix dans shift */
 void Polynoms::shift_diagonal( double *Matrix, const int N, double shift)
 {
 	for (int i=0; i<N; i++)
 		Matrix[i*N + i] += shift;
 };
 //qques outils d'algebre des polynomes

 /** \fn remplir_trinome
  * \brief remplie un trinome
  */
 void Polynoms::remplir_trinome(double coeff, double u, double alpha, vector< double >& trinome)
 {
 	trinome[2] = coeff;
 	coeff *= u;
 	trinome[1] = -2*coeff;
 	coeff *= u;
 	trinome[0] = coeff - alpha;
 };


 /** \fn additionne
  * \brief Additionne 2 vecteurs de tailles differentes
  * @param P,Q les 2 vecteurs à additionner
  * @param n,m les 2 tailles respectives des vecteurs P et Q
  * @return un vecteur qui correspond a P+Q*/
 vector< double > Polynoms::additionne(const vector< double >& P, const vector< double >& Q, int n, int m)
 {
 	if(m < n) {
 		vector< double > P_plus_Q = P;
 		for(int j = 0; j<m+1; j++){
 			P_plus_Q[j] += Q[j];
 		}
 		return P_plus_Q;
 	}
 	else{
 		vector< double > P_plus_Q = Q;
 		for(int j = 0; j<n+1; j++){
 			P_plus_Q[j] += P[j];
 		}
 		return P_plus_Q;
 	}
 };


 /** \fn multiplie
  * \brief Calcule le produit scalaire de 2 vecteurs de tailles differentes
  * @param P,Q les 2 vecteurs
  * @param n,m les 2 tailles respectives des 2 vecteurs
  * @return un vecteur correspond à (P,Q) */
 vector< double > Polynoms::multiplie(const vector< double >& P, const vector< double >& Q, int n, int m)
 {
 	vector< double > P_fois_Q (n+m+1,0.);
 	for(int i = 0; i<n+1; i++){
 		for(int j = 0; j<m+1; j++){
 			P_fois_Q[i+j] += P[i]*Q[j];
 		}
 	}
 	return P_fois_Q;
 }


 /** \fn polynome_caracteristique
  * \brief Calcule le polynome caracteristique
  * @return */
 vector< double > Polynoms::polynome_caracteristique(double alpha_v, double alpha_l, double Uv, double Ul, double rhov, double rhol, double   cv_2, double cl_2, double dPiv, double dPil)//c= inverse vitesse du son!
 {
 	vector< double > trinome_u(3);
 	vector< double > trinome_c(3);
 	vector< double > pol_char;

 	remplir_trinome(alpha_v*cv_2,Uv,alpha_v,trinome_c);
 	remplir_trinome(rhol,Ul,dPil,trinome_u);

 	pol_char = multiplie(trinome_u, trinome_c, 2, 2);

 	remplir_trinome(alpha_l*cl_2,Ul,alpha_l,trinome_c);
 	remplir_trinome(rhov,Uv,dPiv,trinome_u);

 	pol_char = additionne(pol_char,multiplie(trinome_u, trinome_c, 2, 2),4,4);

 	return pol_char;
 }

 /** \fn polynome_caracteristique
  * \brief
  * @param alpha1
  * @param alpha2
  * @param u1
  * @param u2
  * @param rho1
  * @param rho2
  * @param invc1sq
  * @param invc2sq
  * @param dpi1
  * @param dpi2
  * @param g2press
  * @param g2alpha
  * @param g2
  * @param epsilon
  * @return */
 vector< double > Polynoms::polynome_caracteristique(double alpha1, double alpha2, double u1, double u2, double rho1, double rho2, double invc1sq, double invc2sq, double dpi1, double dpi2, double g2press, double g2alpha, double g2, double epsilon)
 {
 	vector< double > trinome1(3);
 	vector< double > trinome2(3);
 	vector< double > pol_char;
 	double K1,K2,K3;

 	/* (x-u1)^2(x-u2)^2-K1(x-u2)^2-K2(x-u1)^2+K3 */
 	K1=alpha1*rho2*g2press + dpi1*alpha2*invc2sq*g2alpha;
 	K2=alpha2*rho1*g2press + dpi2*alpha1*invc1sq*g2alpha;
 	K3=(alpha1*dpi2+alpha2*dpi1)*g2alpha*g2press/g2;
 	//cout<<"K1= " <<K1<<", K2= " <<K2<<", K3= " <<K3<<endl;

 	if(fabs(K1)>epsilon && fabs(K2)>epsilon )
 	{
 		remplir_trinome(1/K1,u1,1,trinome1);
 		//cout<<"coeff constant trinome1= " << trinome1[0]<<endl;
 		remplir_trinome(1/K2,u2,1,trinome2);
 		//cout<<"coeff constant trinome2= " << trinome2[0]<<endl;

 		pol_char = multiplie(trinome1, trinome2, 2, 2);
 		//cout<<"coeff constant produit= " << pol_char[0]<<endl;

 		pol_char[0]+=K3/K2/K1-1;
 		//cout<<"coeff constant apres soustraction= " << pol_char[0]<<endl;
 	}
 	else
 	{
 		remplir_trinome(1,u1,K1,trinome1);
 		remplir_trinome(1,u2,K2,trinome2);

 		pol_char = multiplie(trinome1, trinome2, 2, 2);

 		pol_char[0]+=K3-K2*K1;
 	}

 	/*   for(int ct =0; ct<pol_char.size()-1; ct++)  */
 	/*     if(fabs(pol_char[ct])< epsilon)  */
 	/*       pol_char[ct]=0.; */

 	return pol_char;
 }

 //Pour le tri des valeurs propres

 /** \fn module
  * \brief calcule le module d'un nombre complexe
  * @param z est un nombre complexe
  * @return  calcule le module de z*/
 double Polynoms::module(complex< double > z)
 {
 	return abs(z);
 }

 /** \fn modulec calcule le module² de z */
 double Polynoms::modulec(complex< double > z)
 {
 	return norm(z);
 }

 /** \fn abs_generalise
  * \brief calcule la valeur absolue d'un nombre complexe
  * \Details calcule la valeur absolue d'un nombre complexe en prenant en compte
  * que la partie réelle c-a-d si z= a+ib abs_generalize() renvoie |a|+ib
  * @param z
  * @return si z = a+ib elle renvoie |a|+ib */
 complex< double > Polynoms::abs_generalise(complex< double > z)
 {
 	if(z.real()>=0)
 		return z ;
 	else
 		return -z ;
 }
 void Polynoms::ordre_croissant_abs(vector< complex< double > > &L, int n)
 {
 	vector< complex< double > > copieL = L;
 	int i_max;
 	double Lmaxc;
 	double Lmax_iterc;

 	for(int j=0; j<n; j++){
 		i_max = 0;
 		Lmaxc = copieL[0].real();
 		for (int i=1; i<n; i++){
 			Lmax_iterc = copieL[i].real();
 			if( Lmax_iterc > Lmaxc ){
 				Lmaxc = Lmax_iterc;
 				i_max = i;
 			}
 		}
 		L[n-1-j] = copieL[i_max];
 		copieL[i_max] = -INFINITY;
 	}
 }

 int Polynoms::new_tri_selectif(vector< complex< double > > &L, int n, double epsilon)
 {
 	int size=1;
 	/*
 	cout<<"avant ordre croissant"<<endl;
 	for(int i=0;i<n;i++)
 		cout<<L[i]<<", "<<endl;
 	 */
 	ordre_croissant_abs(L,n);
 	/*
 	cout<<"après ordre croissant"<<endl;
 	for(int i=0;i<n;i++)
 		cout<<L[i]<<", "<<endl;
 	 */
 	vector< complex< double > >result(1,L[0]);
 	for(int i=1;i<n;i++)
 		if(fabs(L[i].real()-result[size-1].real())>epsilon || fabs(L[i].imag())>epsilon)
 		{
 			result.push_back(L[i]);
 			size++;
 			/*
 			cout<<"tri selectif step "<< i<<endl;
 			for(int j=0;j<size;j++)
 				cout<<result[j]<<", "<<endl;
 			 */
 		}
 	L=result;
 	return size;
 }

 /** \fn tri_selectif
  * \brief
  * @param L
  * @param n
  * @param epsilon
  * @return */
 template< typename T >
 int Polynoms::tri_selectif(T& L, int n, double epsilon)
 {
 	int imin;
 	int size_vp =n;
 	int i=0;
 	int j;
 	complex< double > Lmin;

 	for(i = 0 ; i<size_vp ; i++){//On travaille sur le sous vecteur des indices de  i � size_vp-1
 		imin = i;
 		Lmin = L[i];

 		for(j = i+1 ; j < size_vp ; j++)
 			if( Lmin.real() > L[j].real() )
 			{Lmin = L[j]; imin = j;}

 		//on met la vp � zero si elle est trop petite
 		if( fabs(Lmin.real()) < epsilon) Lmin.real()=0;
 		if( fabs(Lmin.imag()) < epsilon) Lmin.imag()=0;

 		switch (i) {
 		case 0: {
 			if (imin != i) {
 				L[imin]=L[i];
 				L[i]=Lmin;
 			}
 			break;
 		}
 		case 1: {
 			if ( module(L[i-1] - Lmin) > epsilon )  {
 				if (imin != i) {
 					L[imin]=L[i];
 					L[i]=Lmin;
 				}

 			} else {
 				L[imin]=L[size_vp - 1];
 				L[size_vp - 1]=Lmin;

 				size_vp--;i--;
 			}
 			break;
 		}
 		default : {
 			if ( (module(L[i-1] - Lmin) > epsilon) &&  (module(L[i-2] - Lmin) > epsilon) )  {
 				if (imin != i) {
 					L[imin]=L[i];
 					L[i]=Lmin;
 				}

 			} else {
 				L[imin]=L[size_vp - 1];
 				L[size_vp - 1]=Lmin;

 				size_vp--;i--;
 			}
 		}
 		}
 	}

 	return size_vp;
 }

 //Calcul des coefficients du polynome d'interpolation x->y dans la base de Newton
 template<class T>

 /** \fn dif_div
  * \brief
  * @param n
  * @param x
  * @param y
  * @param epsilon
  * @return */
 vector< complex< double > > Polynoms::dif_div(int n, const vector< complex< double > >& x, const vector< T >& y, double epsilon)//n=nb valeurs à interpoler
 {
 	complex< double > tab_dif_div[n*n];//La matrice de taille n*n
 	vector< complex< double > > dif_div (n);
 	int i, j;
 	complex< double > delta_x;

 	for( i=0 ; i<n ; i++) {
 		tab_dif_div[i*n]=y[i];
 		//   cerr << tab_dif_div[i*n+0] << endl;
 	}
 	for(j=1; j<n ; j++){
 		for(i=0;i < n-j;i++) {
 			delta_x = x[i+j]-x[i];
 			if(module(delta_x) < epsilon) {
 				cout << " attention, dif_div, delta_x <epsilon, " << " delta_x= " << delta_x <<" epsilon= "<<epsilon<<endl;
 				cout<<" vp= " <<endl;
 				for(int k=0; k<n; k++)
 					cout << x[k]<<", " << endl;
 			}
 			tab_dif_div[i*n+j]=(tab_dif_div[(i+1)*n+j-1]-tab_dif_div[i*n+j-1])/delta_x;
 		}
 	}
 	for(i=0 ; i<n ; i++)
 		dif_div[i] = tab_dif_div[i*n+n-1-i];
 	return dif_div;
 }
 //attention, Les vp complexes conjugu�es doivent se suivre dans la liste x
 void Polynoms::appliquer_dif_div(int n, const vector< complex< double > >& dif_div, const vector< complex< double > >& x, const double* A, const int sizeA, double epsilon, double *p)//p=result
 {
 	int i, debut = 1;
 	double Id[sizeA*sizeA], aux[sizeA*sizeA], matProd[sizeA*sizeA];

 	for(int k=0; k<sizeA*sizeA; k++)
 		Id[k]=0;
 	for(int k=0; k<sizeA; k++)
 		Id[k*sizeA+k]=1;
 	/*   for(int k=0; k<sizeA*sizeA; k++) */
 	/*     aux[k] = A[k] ; */
 	//cerr << x << endl;
 	if( fabs(x[0].imag()) > epsilon){
 		for(int k=0; k<sizeA*sizeA; k++)
 			p[k] = A[k] * dif_div[0].real() ;
 		shift_diagonal(p, sizeA,(dif_div[1] - dif_div[0]*x[1]).real());
 		debut = 2;
 	}else{
 		for(int k=0; k<sizeA*sizeA; k++)
 			p[k] = Id[k] * dif_div[0].real() ;
 	}
 	for(i=debut ; i<n ; i++){
		for(int k=0; k<sizeA*sizeA; k++)
			aux[k] = A[k];
 		//cerr << " on traite la racine reele " << x[i] << endl;
 		if ( fabs(x[i].imag()) < epsilon ) {
 			shift_diagonal(aux, sizeA, (- x[i].real()));
 			matrixProduct( p, sizeA, sizeA, aux, sizeA, sizeA, matProd);
 			for(int k=0; k<sizeA*sizeA; k++)
 				p[k] = matProd[k] ;
 			shift_diagonal(p, sizeA, dif_div[i].real());
 		}
 		else{ //couple de conjugees
 			if(fabs(x[i].imag() + x[i+1].imag()) < epsilon){
 				matrixProduct( aux, sizeA, sizeA, aux, sizeA, sizeA, matProd);
 				for(int k=0; k<sizeA*sizeA; k++)
 					aux[k] = matProd[k];
 				for(int k=0; k<sizeA*sizeA; k++)
 					aux[k]+=(-2*x[i].real()) * A[k];
 				shift_diagonal(aux, sizeA, modulec(x[i]));

 				matrixProduct( p, sizeA, sizeA, aux, sizeA, sizeA, matProd);
 				for(int k=0; k<sizeA*sizeA; k++)
 					p[k] = matProd[k] ;
 				for(int k=0; k<sizeA*sizeA; k++)
 					p[k]+= dif_div[i].real() * A[k];
 				shift_diagonal(p, sizeA, (dif_div[i+1] - dif_div[i]*x[i+1]).real());
 				i++;
 			}
 			else{
 				cout << " appliquer_dif_div : les racines complexes conjuguees ne se suivent pas " <<  endl;
 				for(int k=0; k<n; k++)
 					cout << x[k]<<", " << endl;
 				throw(" appliquer_dif_div : les racines complexes conjuguees ne se suivent pas ");
 			}
 		}
 	}
 }

 void Polynoms::abs_par_interp_directe(int n,  vector< complex< double > > vp,   double * A, int sizeA, double epsilon, double *result,vector< complex< double > > y)
 {
 	vector< complex< double > > differences_divisees = dif_div(n, vp, y, epsilon);
 	/*
 	cout<<"valeurs propres"<<endl;
 	for( int i=0 ; i<n ; i++)
 		cout<<vp[i] <<", "<<endl;
 	cout<<"valeurs à atteindre"<<endl;
 	for( int i=0 ; i<n ; i++)
 		cout<<y[i] <<", "<<endl;
 	cout<<"differences_divisees= "<<endl;
 	for(int i=0 ; i<n ; i++)
 		cout<<differences_divisees[i]<<", ";
 	cout<<endl;
 	 */
 	appliquer_dif_div(n, differences_divisees, vp, A, sizeA, epsilon, result);
 }

 bool Polynoms::belongTo(complex< double > a , vector< complex <double > > v, double epsilon)
 {
 	bool result = false;
 	int size = v.size();
 	int i = 0;
 	double epsilon2 = epsilon*epsilon;

 	while(i<size && !result)
 	{
 		result = modulec(a-v[i]) < epsilon2;
 		i++;
 	}
 	return result;
 }
 double Polynoms::avr(double a, double b)
 {
	 return 0.5*(a+b);
 }
