#ifndef BINVECMUL_H_
#define BINVECMUL_H_

#include "precomp.h"

#define vADD(C,A,B)    {C = _mm_xor_si128(A,B);}
#define vSUB(C,A,B)    vADD(C,A,B)
#define vMULT(C,A,B)   {C = _mm_clmulepi64_si128(A,B,0);}
#define vSFTLB(C,A,B)   {C= _mm_slli_si128(A,B);}
#define vSFTRB(C,A,B)   {C= _mm_srli_si128(A,B);}
#define vAND(C,A,B)    {C= _mm_and_si128(A,B);}


extern void bincopy(gfe1x *c, gfe1x *a);
extern void gfe1xnSq(gfe1x *c, gfe1x *a, int n);
extern void gfe1xMult(gfe1x *c, gfe1x *a, gfe1x *b);
extern void gfe1xMultConst(gfe1x *c, gfe1x *a, vec b);
extern void gfe1xAdd(gfe1x *c, gfe1x *a, gfe1x *b);
extern void ladderStep(gfe1x *x0, gfe1x *x1, gfe1x *z1, gfe1x *x2, gfe1x *z2);
#define gfe1xSq(x,y) gfe1xnSq(x,y,1)


#endif

