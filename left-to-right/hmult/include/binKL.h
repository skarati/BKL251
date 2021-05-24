#ifndef BINKL_H_
#define BINKL_H_
#include "binvecmult.h"

inline void clamp(unsigned char n[32]);
inline void invert(gfe1x *op, gfe1x *in);
extern void conditionalSwap(gfe1x *a, gfe1x *b, gfe1x *c, gfe1x *d, vec sbit);
inline void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]);
inline void scalar_mult_var_base(unsigned char op[32], unsigned char varbase[32], unsigned char n[32]);

inline void clamp(unsigned char n[32]){
	n[0] 	= n[0] & 0xfc;
	n[31] 	= n[31] & 0x07;
	n[31] 	= n[31] | 0x04;
}


void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]){
	int 	prevbit, bit, i, j, k;
	vec 	vbit, t1;
	gfe1x 	sx,sz,rx,rz;
	gfe1x 	t2,t3,t4,t22;
	gfe 	op_gfe;

	
	sx.v[0] = base[0]; sx.v[1]=zero; sx.v[2]=zero; sx.v[3]=zero;
	sz.v[0] = base[1]; sz.v[1]=zero; sz.v[2]=zero; sz.v[3]=zero;
	rx.v[0] = base[2]; rx.v[1]=zero; rx.v[2]=zero; rx.v[3]=zero;
	rz.v[0] = base[3]; rz.v[1]=zero; rz.v[2]=zero; rz.v[3]=zero;

	j = 2;	
	i = 31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1){
			j = 7;
			i--;
		}
	}


	prevbit = 0;
	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = ((n[i]>>j) & 1); 

			//Conditional Swap
			vbit = _mm_set1_epi32(0-(bit^prevbit));
			conditionalSwap(&sx, &rx, &sz, &rz, vbit);
			ladderStepF(&sx, &sz, &rx, &rz);
			
		
			prevbit = bit;
		}
		
		j=7;
	}


	//Conditional Swap
	vbit = _mm_set1_epi32(0-bit);
	conditionalSwap(&sx, &rx, &sz, &rz, vbit);

	invert(&sz,&sz);
	gfe1xMult(&sx, &sx, &sz);

	convert_gfe1x2gfe(&op_gfe, &sx);
	convert_itoc(&op_gfe, op);
}



void scalar_mult_var_base(unsigned char op[32], unsigned char varbase[32], unsigned char n[32]){
	int 	prevbit, bit, i, j, k;
	vec 	vbit, t1;
	gfe1x 	sx,sz,rx,rz,vbx;
	gfe1x 	t2,t3,t4;
	gfe 	op_gfe;
	gfe	work;

	convert_ctoi(&work,varbase);
	convert_gfe2gfe1x(&vbx,&work);

	sx = vbx;
	sz.v[0] = one; sz.v[1]=zero; sz.v[2]=zero; sz.v[3]=zero;
	
	//Set r = KL_Point_Double(s)
	t2.v[0] = _mm_xor_si128(sx.v[0],sz.v[0]);
	t2.v[1]=sx.v[1]; t2.v[2]=sx.v[2]; t2.v[3]=sx.v[3]; //t2 = x1 + z1
	gfe1xnSq(&rx,&t2,2);			//x2 = t2^2
	gfe1xMultConst(&rx, &rx, b);		//x2 = x2*b
	gfe1xSq(&rz,&sx);			//z2 = z2^2
				
	
	j = 2;	
	i = 31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1){
			j = 7;
			i--;
		}
	}


	prevbit = 0;
	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = ((n[i]>>j) & 1); 

			//Conditional Swap
			vbit = _mm_set1_epi32(0-(bit^prevbit));
			conditionalSwap(&sx, &rx, &sz, &rz, vbit);

			ladderStepV(&sx, &sz, &rx, &rz, &vbx);

			prevbit = bit;
		}

		j=7;
	}


	//Conditional Swap
	vbit = _mm_set1_epi32(0-bit);
	conditionalSwap(&sx, &rx, &sz, &rz, vbit);


	invert(&sz,&sz);
	gfe1xMult(&sx, &sx, &sz);

	convert_gfe1x2gfe(&op_gfe, &sx);
	convert_itoc(&op_gfe, op);
}



inline void invert(gfe1x *op, gfe1x *in){
	gfe1x t;
	gfe1x x2,x3,x4,x7,x_6_1,x_12_1,x_24_1,x_25_1,x_50_1;
	gfe1x x_100_1,x_125_1,x_250_1;

	/*2*/			gfe1xSq(&x2, in);				//1S
	/*3*/			gfe1xMult(&x3, &x2, in);			//1M
	/*4*/			gfe1xSq(&x4, &x2);				//2S
	/*7*/			gfe1xMult(&x7, &x4, &x3);			//2M


	/*2^6-3*/		gfe1xnSq(&x_6_1, &x7,3);			//5S
	/*2^6-1*/		gfe1xMult(&x_6_1, &x_6_1,&x7);			//3M
	
	/*2^12-6*/		gfe1xnSq(&x_12_1, &x_6_1,6);			//11S
	/*2^12-1*/		gfe1xMult(&x_12_1, &x_12_1,&x_6_1);		//4M

	/*2^24-1*/		gfe1xnSq(&x_24_1, &x_12_1,12);			//23S	
	/*2^24-1*/		gfe1xMult(&x_24_1, &x_24_1,&x_12_1);		//5M
	
	/*2^25-1*/		gfe1xSq(&x_25_1, &x_24_1);			//24S	
	/*2^25-1*/		gfe1xMult(&x_25_1, &x_25_1,in);			//6M
	
	/*2^50-1*/		gfe1xnSq(&x_50_1, &x_25_1,25);			//49S
	/*2^50-1*/		gfe1xMult(&x_50_1,&x_50_1, &x_25_1);		//7M
	
	/*2^100-1*/		gfe1xnSq(&x_100_1, &x_50_1,50);			//99S
	/*2^100-1*/		gfe1xMult(&x_100_1,&x_100_1, &x_50_1);		//8M
	
	/*2^125-1*/		gfe1xnSq(&x_125_1, &x_100_1,25);		//124S
	/*2^125-1*/		gfe1xMult(&x_125_1,&x_125_1, &x_25_1);		//9M
	
	/*2^250-1*/		gfe1xnSq(&x_250_1, &x_125_1,125);		//249S
	/*2^250-1*/		gfe1xMult(&x_250_1,&x_250_1, &x_125_1);		//10M

	/*2^251-2*/		gfe1xSq(op,&x_250_1);				//250S

}




#endif
