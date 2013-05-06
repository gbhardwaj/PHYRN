/********************************************************/
/*		phyrn_phylip.h				*/
/*******************************************************/
#include <stdio.h>
#include <stdlib.h>

typedef long longer[6];

void initseed(long *inseed, long *inseed0, longer seed)
{ /* input random number seed */
  long i, loopcount;

  assert(inseed);
  assert(inseed0);
  loopcount = 0;

  *inseed0 = *inseed;
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = *inseed & 63;
    *inseed /= 64;
    i++;
  } while (*inseed != 0);
}  /*initseed*/

double randum(longer seed)
{ /* random number generator -- slow but machine independent
  This is a multiplicative congruential 32-bit generator
  x(t+1) = 1664525 * x(t) mod 2^32, one that passes the
  Coveyou-Macpherson and Lehmer tests, see Knuth ACP vol. 2
  We here implement it representing each integer in base-64
  notation -- i.e. as an array of 6 six-bit chunks   */

  long i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;   /* these four statements set the multiplier */
  mult[1] = 24;   /* -- they are its "digits" in a base-64    */
  mult[2] = 22;   /*    notation: 1664525 = 13*64^3+24*64^2   */
  mult[3] = 6;    /*                         +22*64+6         */
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {  /* do the multiplication piecewise */
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));        /* new seed replaces old one */
  seed[5] &= 3;          /* from the new seed, get a floating point fraction */
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */

