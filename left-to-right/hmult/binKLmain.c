#include "basic.h"
#include "binvecmult.h"
#include "binKL.h"
#include "measurement.h"

int main(){
	int de;
	char vbase[32];
	char np[32];
	char npp[32];
	


	char n[32] = {183, 201,  60, 253, 218, 90,  195, 130, 
                      133, 189, 114,  23,  21, 253,  99,  63, 
                      163, 204, 217,  98, 161,  83,  41, 110, 
                      242, 191,  10, 193, 252,  61,  49, 1};

	clamp(n);


	MEASURE({
		scalar_mult_fixed_base(np, n);
	});
	printf("Total CPU cycles for fixed-base scalar multiplication Min: %.2f.\n", RDTSC_clk_min);
	printf("Total CPU cycles for fixed-base scalar multiplication Median: %.2f.\n", RDTSC_clk_median);
	printf("Total CPU cycles for fixed-base scalar multiplication Max: %.2f.\n", RDTSC_clk_max);
	printf("\n\n");

	de = 255;
	printf("\nnp= ");
	for(int i= 31; i>=0; i--){
		//printf("\ni=%d\n",i);
		for(int j=7;j>=0;j--){
			if(((np[i]>>j)&1)==1)
				printf("z^%d +",de);
			de--;
		}
	}	

	printf("\n\n\n");
	MEASURE({
		scalar_mult_var_base(npp, np, n);
	});
	printf("Total CPU cycles for variable-base scalar multiplication Min: %.2f.\n", RDTSC_clk_min);
	printf("Total CPU cycles for variable-base scalar multiplication Median: %.2f.\n", RDTSC_clk_median);
	printf("Total CPU cycles for variable-base scalar multiplication Max: %.2f.\n", RDTSC_clk_max);
	printf("\n\n");

	de = 255;
	printf("\nnpp= ");
	for(int i= 31; i>=0; i--){
		//printf("\ni=%d\n",i);
		for(int j=7;j>=0;j--){
			if(((npp[i]>>j)&1)==1)
				printf("z^%d +",de);
			de--;
		}
	}	

	printf("\n\n\n");




	return 0;
}

