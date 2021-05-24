#ifndef BASICS_H_
#define BASICS_H_

#include <immintrin.h>
#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <inttypes.h>
#include <stdint.h>


#define align __attribute__ ((aligned (32)))

#if !defined (ALIGN16)
	# if defined (__GNUC__)
		# define ALIGN32 __attribute__ ( (aligned (32)))
	# else
		# define ALIGN32 __declspec (align (32))
	# endif
#endif

#define CRYPTO_SCALARBYTES 32
#define CRYPTO_BYTES 48

typedef unsigned char u8;
typedef unsigned int u32;
typedef unsigned long long u64;
typedef unsigned __int128 u128;
typedef __m128i vec;

typedef struct{
	u64 v[4];
}gfe;

typedef struct{
	u64 v[5];
}gfe51;

typedef struct {
  	vec v[4];
}gfe1x;

const vec vMask64 = {0xffffffffffffffff,0x0};
const vec vMask59 = {0x7ffffffffffffff,0x0};
const vec zero = {0x0,0x0};
const vec one = {0x1,0x0};
const vec zeroOne = {-1,0x0};
const vec residue1 = {0x95,0x0};
const vec residue2 = {0x12A0,0x0};
const vec b = {0x2387,0x0};
const vec baseFixed = {0xC,0x0};
const vec base[4] = {{0xC,0x0},{0x1,0x0},{0x21BD487,0x0},{0x50,0x0}};

void convert_ctoi(gfe *r64, const unsigned char r[32]){
	r64->v[0]  = r[0];
	r64->v[0] |= ((u64)r[1] << 8);
	r64->v[0] |= ((u64)r[2] << 16);
	r64->v[0] |= ((u64)r[3] << 24);
	r64->v[0] |= ((u64)r[4] << 32);
	r64->v[0] |= ((u64)r[5] << 40);
	r64->v[0] |= ((u64)r[6] << 48);
	r64->v[0] |= ((u64)r[7] << 56);

	r64->v[1]  = r[8];
	r64->v[1] |= ((u64)r[9] << 8);
	r64->v[1] |= ((u64)r[10] << 16);
	r64->v[1] |= ((u64)r[11] << 24);
	r64->v[1] |= ((u64)r[12] << 32);
	r64->v[1] |= ((u64)r[13] << 40);
	r64->v[1] |= ((u64)r[14] << 48);
	r64->v[1] |= ((u64)r[15] << 56);

	r64->v[2]  = r[16];
	r64->v[2] |= ((u64)r[17] << 8);
	r64->v[2] |= ((u64)r[18] << 16);
	r64->v[2] |= ((u64)r[19] << 24);
	r64->v[2] |= ((u64)r[20] << 32);
	r64->v[2] |= ((u64)r[21] << 40);
	r64->v[2] |= ((u64)r[22] << 48);
	r64->v[2] |= ((u64)r[23] << 56);

	r64->v[3]  = r[24];
	r64->v[3] |= ((u64)r[25] << 8);
	r64->v[3] |= ((u64)r[26] << 16);
	r64->v[3] |= ((u64)r[27] << 24);
	r64->v[3] |= ((u64)r[28] << 32);
	r64->v[3] |= ((u64)r[29] << 40);
	r64->v[3] |= ((u64)r[30] << 48);
	r64->v[3] |= ((u64)(r[31]&0x7) << 56);
}




void convert_itoc(gfe *r64, unsigned char r[32]){
	r[0] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[1] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[2] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[3] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[4] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[5] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[6] = r64->v[0] & 0xff;
	r64->v[0] = r64->v[0]>>8;
	r[7] = r64->v[0] & 0xff;

	r[8] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[9] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[10] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[11] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[12] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[13] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[14] = r64->v[1] & 0xff;
	r64->v[1] = r64->v[1]>>8;
	r[15] = r64->v[1] & 0xff;

	r[16] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[17] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[18] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[19] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[20] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[21] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[22] = r64->v[2] & 0xff;
	r64->v[2] = r64->v[2]>>8;
	r[23] = r64->v[2] & 0xff;

	r[24] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[25] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[26] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[27] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[28] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[29] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[30] = r64->v[3] & 0xff;
	r64->v[3] = r64->v[3]>>8;
	r[31] = r64->v[3] & 0x7;
}

void convert_gfe1x2gfe(gfe *op, gfe1x *in){
	op->v[0] = 0[(u64 *) &in->v[0]];
	op->v[1] = 0[(u64 *) &in->v[1]];
	op->v[2] = 0[(u64 *) &in->v[2]];
	op->v[3] = 0[(u64 *) &in->v[3]];
	
}

void convert_gfe2gfe1x(gfe1x *op, gfe *in){
	 0[(u64 *) &op->v[0]]= in->v[0];	1[(u64 *) &op->v[0]]= 0;
	 0[(u64 *) &op->v[1]]= in->v[1];	1[(u64 *) &op->v[1]]= 0;
	 0[(u64 *) &op->v[2]]= in->v[2];	1[(u64 *) &op->v[2]]= 0;
	 0[(u64 *) &op->v[3]]= in->v[3];	1[(u64 *) &op->v[3]]= 0;
	
}


#endif


