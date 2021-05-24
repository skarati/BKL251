.p2align 5
.globl gfe1xMult
gfe1xMult:

	movdqa	0(%rsi),%xmm8			#store first 16 bytes from second argument to xmm8
	movdqa	16(%rsi),%xmm9			#store second 16 bytes from second argument to xmm9
	movdqa	32(%rsi),%xmm10			#store third 16 bytes from second argument to xmm10
	movdqa	48(%rsi),%xmm11			#store forth 16 bytes from second argument to xmm11

	movdqa	0(%rdx),%xmm12			#store first 16 bytes from third argument to xmm12
	movdqa	16(%rdx),%xmm13			#store second 16 bytes from third argument to xmm13
	movdqa	32(%rdx),%xmm14			#store third 16 bytes from third argument to xmm14
	movdqa	48(%rdx),%xmm15			#store forth 16 bytes from third argument to xmm15

	movdqa	%xmm12, %xmm0			#xmm0 = xmm12
	pclmulqdq $0x00, %xmm8, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm13, %xmm1			#xmm1 = xmm13
	pclmulqdq $0x00, %xmm8, %xmm1		#xmm1 = xmm8 * xmm13
	movdqa	%xmm14, %xmm2			#xmm2 = xmm14
	pclmulqdq $0x00, %xmm8, %xmm2		#xmm2 = xmm8 * xmm14
	movdqa	%xmm15, %xmm3			#xmm3 = xmm15
	pclmulqdq $0x00, %xmm8, %xmm3		#xmm3 = xmm8 * xmm15

	movdqa	%xmm12, %xmm7			#xmm7 = xmm12
	pclmulqdq $0x00, %xmm9, %xmm7		#xmm7 = xmm9 * xmm12
	pxor %xmm7, %xmm1			#xmm1 = xmm8 * xmm13 + xmm9 * xmm12
	movdqa	%xmm13, %xmm7			#xmm7 = xmm13
	pclmulqdq $0x00, %xmm9, %xmm7		#xmm7 = xmm9 * xmm13
	pxor %xmm7, %xmm2			#xmm2 = xmm8 * xmm14 + xmm9 * xmm13
	movdqa	%xmm14, %xmm7			#xmm7 = xmm14
	pclmulqdq $0x00, %xmm9, %xmm7		#xmm7 = xmm9 * xmm14
	pxor %xmm7, %xmm3			#xmm3 = xmm8 * xmm15 + xmm9 * xmm14
	movdqa	%xmm15, %xmm4			#xmm4 = xmm15
	pclmulqdq $0x00, %xmm9, %xmm4		#xmm4 = xmm9 * xmm15


	movdqa	%xmm12, %xmm7			#xmm7 = xmm12
	pclmulqdq $0x00, %xmm10, %xmm7		#xmm7 = xmm10 * xmm12
	pxor %xmm7, %xmm2			#xmm2 = xmm8 * xmm14 + xmm9 * xmm13 + xmm10 * xmm12
	movdqa	%xmm13, %xmm7			#xmm7 = xmm13
	pclmulqdq $0x00, %xmm10, %xmm7		#xmm7 = xmm10 * xmm13
	pxor %xmm7, %xmm3			#xmm3 = xmm8 * xmm15 + xmm9 * xmm14 + xmm10 * xmm13
	movdqa	%xmm14, %xmm7			#xmm7 = xmm14
	pclmulqdq $0x00, %xmm10, %xmm7		#xmm7 = xmm10 * xmm14
	pxor %xmm7, %xmm4			#xmm4 = xmm9 * xmm15 + xmm10 * xmm14
	movdqa	%xmm15, %xmm5			#xmm5 = xmm15
	pclmulqdq $0x00, %xmm10, %xmm5		#xmm5 = xmm10 * xmm15

	movdqa	%xmm12, %xmm7			#xmm7 = xmm12
	pclmulqdq $0x00, %xmm11, %xmm7		#xmm7 = xmm11 * xmm12
	pxor %xmm7, %xmm3			#xmm3 = xmm8 * xmm15 + xmm9 * xmm14 + xmm10 * xmm13 + xmm11 * xmm12
	movdqa	%xmm13, %xmm7			#xmm7 = xmm13
	pclmulqdq $0x00, %xmm11, %xmm7		#xmm7 = xmm11 * xmm13
	pxor %xmm7, %xmm4			#xmm4 = xmm9 * xmm15 + xmm10 * xmm14 + xmm11 * xmm13
	movdqa	%xmm14, %xmm7			#xmm7 = xmm14
	pclmulqdq $0x00, %xmm11, %xmm7		#xmm7 = xmm11 * xmm14
	pxor %xmm7, %xmm5			#xmm5 = xmm10 * xmm15 + xmm11 * xmm14
	movdqa	%xmm15, %xmm6			#xmm6 = xmm15
	pclmulqdq $0x00, %xmm11, %xmm6		#xmm6 = xmm11 * xmm15

	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,0(%rdi)			#store xmm0 to first 16 bytes from first argument
	movdqa	%xmm1,16(%rdi)			#store xmm1 to second 16 bytes from first argument
	movdqa	%xmm2,32(%rdi)			#store xmm2 to third 16 bytes from first argument
	movdqa	%xmm3,48(%rdi)			#store xmm3 to forth 16 bytes from first argument

ret
