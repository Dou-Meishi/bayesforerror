/* filename: _func.h.  */

#define SIGMA 1.0
#define CBAR_LE .001
#define CBAR_GE 1000
#define RAND_COUNT 1000


/* function prototypes */

double _theta(double);

double _pr_cbar_A(double);
double _pr_cbar_B(double);
double _pr_cbar_C(double);

double _pr_cn_if_cbar_A(double cn, double cbar);
double _pr_cn_if_cbar_B(double cn, double cbar);
double _pr_cn_if_cbar_C(double cn, double cbar);

double _A_delta_if_cbar_f(double t, double delta,  double cbar,
                                 double Q, int h, int k);

double _borwein_volQ(double a0, const double *ak, int n);
double _fct_sinc(double delta, double cbar, double Q, int h, int k);
double _pr_h_delta_cbar(double delta, double cbar, double Q, int h, int k);
