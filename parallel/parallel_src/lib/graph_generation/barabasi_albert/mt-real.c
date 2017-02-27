/* A C-program for MT19937B: Real number version */
/* genrand() generate one pseudorandom number with double precision   */
/* which is uniformly distributed on [0,1]-interval for each call.    */
/* sgenrand(seed) set initial values to the working area of 624 words.*/
/* sgenrand(seed) must be called once before calling genrand()        */
/* (seed is any integer except 0).                                    */

/* 
    LICENCE CONDITIONS: 

		Matsumoto and Nishimura consent to GNU General 
		Public Licence

    NOTE: 
		When you use it in your program, please let Matsumoto
                <matumoto@math.keio.ac.jp> know it.

		Because of a machine-trouble, Matsumoto lost emails 
   		which arrived during May 28-29.
*/


#include<stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* for tempering */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long ptgfsr[N]; /* set initial seeds: N = 624 words */

void
sgenrand(unsigned long seed) /* seed should not be 0 */
{
  int k;
  
  /* setting initial seeds to ptgfsr[N] using     */
  /* the generator Line 25 of Table 1 in          */
  /* [KNUTH 1981, The Art of Computer Programming */
  /*    Vol. 2 (2nd Ed.), pp102]                  */

  ptgfsr[0]= seed & 0xffffffff;
  for (k=1; k<N; k++)
	ptgfsr[k] = (69069 * ptgfsr[k-1]) & 0xffffffff;
}

double
genrand()
{
  unsigned long y;
  static int k = 1;
  static unsigned long mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if(k == N){ /* generate N words at one time */
	int kk;
	for (kk=0;kk<N-M;kk++) {
	  y = (ptgfsr[kk]&UPPER_MASK)|(ptgfsr[kk+1]&LOWER_MASK);
	  ptgfsr[kk] = ptgfsr[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for (;kk<N-1;kk++) {
	  y = (ptgfsr[kk]&UPPER_MASK)|(ptgfsr[kk+1]&LOWER_MASK);
	  ptgfsr[kk] = ptgfsr[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (ptgfsr[N-1]&UPPER_MASK)|(ptgfsr[0]&LOWER_MASK);
	ptgfsr[N-1] = ptgfsr[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
	
	k = 0;
  }
  
  y = ptgfsr[k++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y &= 0xffffffff; /* you may delete this line if word size = 32 */
  y ^= TEMPERING_SHIFT_L(y);

  return ( (double)y / (unsigned long)0xffffffff );
}


