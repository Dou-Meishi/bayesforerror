/* filename: _func.c */

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "_funcs.h"


/* function implementations */

static const double _1_on_sqrt_2pi = .39894228040143267794;
static const double _1_on_log_gele = .07238241365054197127;


double _theta(double x){
  return (x > 0) ? 1 : 0;
}

double _pr_cbar_A(double cbar){
  return _1_on_log_gele/cbar * _theta(cbar-CBAR_LE)*_theta(CBAR_GE-cbar);
}

double _pr_cbar_B(double cbar){
  return _1_on_sqrt_2pi/(cbar*SIGMA) *
    exp(-log(cbar)*log(cbar)*.5/(SIGMA*SIGMA));
}

double _pr_cbar_C(double cbar){
  return _1_on_log_gele/cbar * _theta(cbar-CBAR_LE)*_theta(CBAR_GE-cbar);
}

double _pr_cn_if_cbar_A(double cn, double cbar){
  return .5/cbar * _theta(cbar - fabs(cn));
}

double _pr_cn_if_cbar_B(double cn, double cbar){
  return .5/cbar * _theta(cbar - fabs(cn));
}

double _pr_cn_if_cbar_C(double cn, double cbar){
  return _1_on_sqrt_2pi/cbar * exp(-cn*cn*.5/(cbar*cbar));
}

double _A_delta_if_cbar_f(double t, double delta, double cbar,
                       double Q, int h, int k)
{
  double res = cos(delta * t);
  double ait = pow(Q, k+1) * t * cbar;

  /*
  printf("%f\n", ait);
  printf("t: %f, delta: %f, cbar: %f, Q: %f, h: %d k: %d\n",
  t, delta, cbar, Q, h, k); */

  for (int i=0; i < h; i++)
    {
      if ( fabs(ait) > 1e-8 ) {
        res *= sin(ait) / ait;
        ait *= Q;
      }
      else
        break;
    }

  return res;
}


/* return the integral of sinc a0 * prod(sinc ak) on (-inf,+inf),
   unit is pi/a0.  */
double _borwein_volQ(double a0, const double *ak, int n)
{
  if (n==0)
    return a0;                   /* probability wont greater than 1 */
  else
    ;
  
  int count = 0;
  double sum_res = 0;
  
  for (int i=0; i < RAND_COUNT; i++)
    {
      sum_res = 0;
      for (int k=0; k < n; k++)
        sum_res += ((double) rand()/RAND_MAX *2 - 1) * ak[k];

      if ( fabs(sum_res) < a0)
        count ++;
      else
        ;
    }

  return (double) count/RAND_COUNT;
}

/* return the integral of PLB2019 Eq(7). */
double _fct_sinc(double delta, double cbar, double Q, int h, int k)
{
  double res = 0;
  double _a0 = cbar * pow(Q, k+1);
  double *am = (double *) malloc((h-1) * sizeof(double)); /* h>=1 */

  if (am == NULL)
    ;                           /* case of h=1 */
  else {
    am[0] = _a0*Q;
    for (int m=1; m < h-1; m++)
      am[m] = am[m-1] * Q;
    }

  res += _borwein_volQ(_a0+delta, am, h-1);
  
  if ( _a0 - delta > 0)
    res += _borwein_volQ(_a0-delta, am, h-1);
  else
    res -= _borwein_volQ(delta-_a0, am, h-1);

  if (am == NULL)
    ;
  else
    free(am);
  
  return res*.25/_a0;
}
