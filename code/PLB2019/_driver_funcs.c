#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "_funcs.h"

int main(void)
{
  printf("%f\n", _A_delta_if_cbar_f(1, 0, 1, 1, 1, 1));
  printf("%f\n", _A_delta_if_cbar_f(1, 0, 1, 1, 2, 1));
  putchar('\n');

  /* time_t t; */
  /* srand((unsigned) time(&t)); */
  /* by default, srand(1) is called. */

  printf("%f\n", _borwein_volQ(1., NULL, 0));
  printf("%f\n", _borwein_volQ(1., (double []) {1.}, 1));
  printf("%f\n", _borwein_volQ(1., (double []) {2.}, 1));
  printf("%f\n", _borwein_volQ(2., (double []) {1.}, 1));
  printf("%f\n", _borwein_volQ(1., (double []) {1/3, 1/5, 1/7}, 3));
  printf("%f\n", _borwein_volQ(2., NULL, 0));
  putchar('\n');

  printf("%f\n", _fct_sinc(0, 1, 1, 1, 1));
  printf("%f\n", _fct_sinc(0, 1, 1, 2, 1));
  printf("%f\n", _fct_sinc(2, 1, 1, 1, 1));

  /* 
  for (int i=0; i < 1000; i++)
    _fct_sinc(.001, .1, .5, 10, 2);
  */

  putchar('\n');
  void test_itos(int);

  for (int num=0; num < 8; num++)
    test_itos(num);

  return 0;
}

void test_itos(int num)
{
  short * s = (short *) malloc(4*sizeof(short));

  if (s != NULL){
    for (int i=3; i>=0; i--, num >>= 1)
      s[i] = ( num & 1 ) ? 1 : -1;
  }
  else
    ;

  for (int i=0; i < 4; i++)
    printf("%3.hd", s[i]);
  putchar('\n');

  free(s);
}
