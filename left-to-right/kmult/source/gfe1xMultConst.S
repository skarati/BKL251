.p2align 5
.globl gfe1xMultConst
gfe1xMultConst:

	movdqa	0(%rsi),%xmm1			#store first 16 bytes from second argument to xmm1
	movdqa	16(%rsi),%xmm2			#store second 16 bytes from second argument to xmm2
	movdqa	32(%rsi),%xmm3			#store third 16 bytes from second argument to xmm3
	movdqa	48(%rsi),%xmm4			#store forth 16 bytes from second argument to xmm4

	#third argument "vec b" is already in xmm0

	pclmulqdq $0x00, %xmm0, %xmm1		#xmm1 = xmm1 * xmm0
	pclmulqdq $0x00, %xmm0, %xmm2		#xmm2 = xmm2 * xmm0
	pclmulqdq $0x00, %xmm0, %xmm3		#xmm3 = xmm3 * xmm0
	pclmulqdq $0x00, %xmm0, %xmm4		#xmm4 = xmm4 * xmm0

	#final reduction
	movdqa %xmm1, %xmm0			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm4 = xmm4 >> 64
	pxor %xmm1, %xmm2			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm2, %xmm1			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm2			#xmm5 = xmm5 >> 64
	pxor %xmm2, %xmm3			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm3, %xmm2			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm6 = xmm6 >> 64
	pxor %xmm3, %xmm4			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm4, %xmm3			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm4			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm4				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm4	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,0(%rdi)			#store xmm0 to first 16 bytes from first argument
	movdqa	%xmm1,16(%rdi)			#store xmm1 to second 16 bytes from first argument
	movdqa	%xmm2,32(%rdi)			#store xmm2 to third 16 bytes from first argument
	movdqa	%xmm3,48(%rdi)			#store xmm3 to forth 16 bytes from first argument

ret
