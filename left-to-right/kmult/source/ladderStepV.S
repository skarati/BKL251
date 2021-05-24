.p2align 5
.globl ladderStepV
ladderStepV:

	movq   %rsp, %r11			#r11 = rsp
	subq   $400, %rsp			#create a stack of 256 byte
	movq   %r11,  0(%rsp)			#push r11 in to stack to restore rsp at the end

	
	###############
	#t22 = (x1+z1)#
	###############
	movdqa	0(%rdi),%xmm4			#store first 16 bytes from first argument to xmm4
	movdqa	16(%rdi),%xmm5			#store second 16 bytes from first argument to xmm5
	movdqa	32(%rdi),%xmm6			#store third 16 bytes from first argument to xmm6
	movdqa	48(%rdi),%xmm7			#store forth 16 bytes from first argument to xmm7

	movdqa	0(%rsi),%xmm8			#store first 16 bytes from second argument to xmm8
	movdqa	16(%rsi),%xmm9			#store second 16 bytes from second argument to xmm9
	movdqa	32(%rsi),%xmm10			#store third 16 bytes from second argument to xmm10
	movdqa	48(%rsi),%xmm11			#store forth 16 bytes from second argument to xmm11

	pxor %xmm4, %xmm8			#xmm8 = xmm8 ^ xmm4 = xmm8 + xmm4
	pxor %xmm5, %xmm9			#xmm9 = xmm9 ^ xmm5 = xmm9 + xmm5
	pxor %xmm6, %xmm10			#xmm10 = xmm10 ^ xmm6 = xmm10 + xmm6
	pxor %xmm7, %xmm11			#xmm11 = xmm11 ^ xmm7 = xmm11 + xmm7

	movdqa	%xmm8,8(%rsp)			#store t22 to stack
	movdqa	%xmm9,24(%rsp)			#store t22 to stack
	movdqa	%xmm10,40(%rsp)			#store t22 to stack
	movdqa	%xmm11,56(%rsp)			#store t22 to stack



	##############
	#t3 = (x2+z2)#
	##############
	movdqa	0(%rdx),%xmm4			#store first 16 bytes from third argument to xmm4
	movdqa	16(%rdx),%xmm5			#store second 16 bytes from third argument to xmm5
	movdqa	32(%rdx),%xmm6			#store third 16 bytes from third argument to xmm6
	movdqa	48(%rdx),%xmm7			#store forth 16 bytes from third argument to xmm7

	movdqa	0(%rcx),%xmm12			#store first 16 bytes from fourth argument to xmm12
	movdqa	16(%rcx),%xmm13			#store second 16 bytes from fourth argument to xmm13
	movdqa	32(%rcx),%xmm14			#store third 16 bytes from fourth argument to xmm14
	movdqa	48(%rcx),%xmm15			#store forth 16 bytes from fourth argument to xmm15

	pxor %xmm4, %xmm12			#xmm12 = xmm12 ^ xmm4 = xmm12 + xmm4
	pxor %xmm5, %xmm13			#xmm13 = xmm13 ^ xmm5 = xmm13 + xmm5
	pxor %xmm6, %xmm14			#xmm14 = xmm14 ^ xmm6 = xmm14 + xmm6
	pxor %xmm7, %xmm15			#xmm15 = xmm15 ^ xmm7 = xmm15 + xmm7


	######################
	#t4 = (x1+z1)*(x2+z2)#
	######################
	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9

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


	###########
	#t4 = t4^2#
	###########
	movdqa %xmm3, %xmm6			#xmm6 = xmm3
	movdqa %xmm2, %xmm4			#xmm4 = xmm2
	movdqa %xmm1, %xmm2			#xmm2 = xmm1

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
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


	movdqa	%xmm0,72(%rsp)			#store t4 to stack
	movdqa	%xmm1,88(%rsp)			#store t4 to stack
	movdqa	%xmm2,104(%rsp)			#store t4 to stack
	movdqa	%xmm3,120(%rsp)			#store t4 to stack


	############
	#t2 = x1*z2#
	############
	movdqa	0(%rdi),%xmm8			#store first 16 bytes from first argument to xmm8
	movdqa	16(%rdi),%xmm9			#store second 16 bytes from first argument to xmm9
	movdqa	32(%rdi),%xmm10			#store third 16 bytes from first argument to xmm10
	movdqa	48(%rdi),%xmm11			#store forth 16 bytes from first argument to xmm11

	movdqa	0(%rcx),%xmm12			#store first 16 bytes from fourth argument to xmm12
	movdqa	16(%rcx),%xmm13			#store second 16 bytes from fourth argument to xmm13
	movdqa	32(%rcx),%xmm14			#store third 16 bytes from fourth argument to xmm14
	movdqa	48(%rcx),%xmm15			#store forth 16 bytes from fourth argument to xmm15

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9

	movdqa	%xmm0,136(%rsp)			#store t2 to stack
	movdqa	%xmm1,152(%rsp)			#store t2 to stack
	movdqa	%xmm2,168(%rsp)			#store t2 to stack
	movdqa	%xmm3,184(%rsp)			#store t2 to stack
	movdqa	%xmm4,200(%rsp)			#store t2 to stack
	movdqa	%xmm5,216(%rsp)			#store t2 to stack
	movdqa	%xmm6,232(%rsp)			#store t2 to stack


	

	
	############
	#t3 = x2*z1#
	############
	movdqa	0(%rsi),%xmm8			#store first 16 bytes from second argument to xmm8
	movdqa	16(%rsi),%xmm9			#store second 16 bytes from second argument to xmm9
	movdqa	32(%rsi),%xmm10			#store third 16 bytes from second argument to xmm10
	movdqa	48(%rsi),%xmm11			#store forth 16 bytes from second argument to xmm11

	movdqa	0(%rdx),%xmm12			#store first 16 bytes from third argument to xmm12
	movdqa	16(%rdx),%xmm13			#store second 16 bytes from third argument to xmm13
	movdqa	32(%rdx),%xmm14			#store third 16 bytes from third argument to xmm14
	movdqa	48(%rdx),%xmm15			#store forth 16 bytes from third argument to xmm15

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9



	####################
	#z2 = x1*z2 + x2*z1#
	####################

	movdqa	136(%rsp),%xmm7			#restore t2 from stack
	movdqa	152(%rsp),%xmm8			#restore t2 from stack
	movdqa	168(%rsp),%xmm9			#restore t2 from stack
	movdqa	184(%rsp),%xmm10		#restore t2 from stack
	movdqa	200(%rsp),%xmm11		#restore t2 from stack
	movdqa	216(%rsp),%xmm12		#restore t2 from stack
	movdqa	232(%rsp),%xmm13		#restore t2 from stack

	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7 = xmm0 + xmm7
	pxor %xmm8, %xmm1			#xmm1 = xmm1 ^ xmm8 = xmm1 + xmm8
	pxor %xmm9, %xmm2			#xmm2 = xmm2 ^ xmm9 = xmm2 + xmm9
	pxor %xmm10, %xmm3			#xmm3 = xmm3 ^ xmm10 = xmm3 + xmm10
	pxor %xmm11, %xmm4			#xmm4 = xmm4 ^ xmm11 = xmm4 + xmm11
	pxor %xmm12, %xmm5			#xmm5 = xmm5 ^ xmm12 = xmm5 + xmm12
	pxor %xmm13, %xmm6			#xmm6 = xmm6 ^ xmm13 = xmm6 + xmm13

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

	###########
	#z2 = z2^2#
	###########
	movdqa %xmm3, %xmm6			#xmm6 = xmm3
	movdqa %xmm2, %xmm4			#xmm4 = xmm2
	movdqa %xmm1, %xmm2			#xmm2 = xmm1

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
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

	
	##############
	#x2 = t4 - z2#
	##############
	movdqa	72(%rsp),%xmm4			#restore t4 from stack
	movdqa	88(%rsp),%xmm5			#restore t4 from stack
	movdqa	104(%rsp),%xmm6			#restore t4 from stack
	movdqa	120(%rsp),%xmm7			#restore t4 from stack

	pxor %xmm0,%xmm4			#xmm4 = xmm0 ^ xmm4 = xmm0 + xmm4
	pxor %xmm1,%xmm5			#xmm5 = xmm1 ^ xmm5 = xmm1 + xmm5
	pxor %xmm2,%xmm6			#xmm6 = xmm2 ^ xmm6 = xmm2 + xmm6
	pxor %xmm3,%xmm7			#xmm7 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa	%xmm4,0(%rdx)			#store xmm0 to first 16 bytes from third argument
	movdqa	%xmm5,16(%rdx)			#store xmm1 to second 16 bytes from third argument
	movdqa	%xmm6,32(%rdx)			#store xmm2 to third 16 bytes from third argument
	movdqa	%xmm7,48(%rdx)			#store xmm3 to forth 16 bytes from third argument


	###########
	#z2 = z2*vbx#
	###########
	movdqa	%xmm0,%xmm8			#xmm8 = xmm0
	movdqa	%xmm1,%xmm9			#xmm9 = xmm1
	movdqa	%xmm2,%xmm10			#xmm10 = xmm2
	movdqa	%xmm3,%xmm11			#xmm11 = xmm3

	movdqa	0(%r8),%xmm12			#store first 16 bytes from third argument to xmm12
	movdqa	16(%r8),%xmm13			#store second 16 bytes from third argument to xmm13
	movdqa	32(%r8),%xmm14			#store third 16 bytes from third argument to xmm14
	movdqa	48(%r8),%xmm15			#store forth 16 bytes from third argument to xmm15

	#movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9

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

	movdqa	%xmm0,0(%rcx)			#store xmm0 to first 16 bytes from fourth argument
	movdqa	%xmm1,16(%rcx)			#store xmm1 to second 16 bytes from fourth argument
	movdqa	%xmm2,32(%rcx)			#store xmm2 to third 16 bytes from fourth argument
	movdqa	%xmm3,48(%rcx)			#store xmm3 to forth 16 bytes from fourth argument

	############
	#z1 = x1*z1#
	############
	movdqa	0(%rdi),%xmm8			#store first 16 bytes from first argument to xmm8
	movdqa	16(%rdi),%xmm9			#store second 16 bytes from first argument to xmm9
	movdqa	32(%rdi),%xmm10			#store third 16 bytes from first argument to xmm10
	movdqa	48(%rdi),%xmm11			#store forth 16 bytes from first argument to xmm11

	movdqa	0(%rsi),%xmm12			#store first 16 bytes from second argument to xmm12
	movdqa	16(%rsi),%xmm13			#store second 16 bytes from second argument to xmm13
	movdqa	32(%rsi),%xmm14			#store third 16 bytes from second argument to xmm14
	movdqa	48(%rsi),%xmm15			#store forth 16 bytes from second argument to xmm15

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9

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


	###########
	#z1 = z1^2#
	###########
	movdqa %xmm3, %xmm6			#xmm6 = xmm3
	movdqa %xmm2, %xmm4			#xmm4 = xmm2
	movdqa %xmm1, %xmm2			#xmm2 = xmm1

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
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

	movdqa	%xmm0,0(%rsi)			#store xmm0 to first 16 bytes from second argument
	movdqa	%xmm1,16(%rsi)			#store xmm1 to second 16 bytes from second argument
	movdqa	%xmm2,32(%rsi)			#store xmm2 to third 16 bytes from second argument
	movdqa	%xmm3,48(%rsi)			#store xmm3 to forth 16 bytes from second argument

	################
	#x1 = (x1+z1)^4#
	################
	movdqa	8(%rsp),%xmm0			#restore (x1+z1) from stack
	movdqa	24(%rsp),%xmm2			#restore (x1+z1) from stack
	movdqa	40(%rsp),%xmm4			#restore (x1+z1) from stack
	movdqa	56(%rsp),%xmm6			#restore (x1+z1) from stack

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
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


	movdqa %xmm3, %xmm6			#xmm6 = xmm3
	movdqa %xmm2, %xmm4			#xmm4 = xmm2
	movdqa %xmm1, %xmm2			#xmm2 = xmm1

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
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

	
	###########
	#x1 = x1*b#
	###########
	pclmulqdq $0x00, b, %xmm0		#xmm0 = xmm0 * b
	pclmulqdq $0x00, b, %xmm1		#xmm1 = xmm1 * b
	pclmulqdq $0x00, b, %xmm2		#xmm2 = xmm2 * b
	pclmulqdq $0x00, b, %xmm3		#xmm3 = xmm3 * b

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

	movdqa	%xmm0,0(%rdi)			#store xmm0 to first 16 bytes from second argument
	movdqa	%xmm1,16(%rdi)			#store xmm1 to second 16 bytes from second argument
	movdqa	%xmm2,32(%rdi)			#store xmm2 to third 16 bytes from second argument
	movdqa	%xmm3,48(%rdi)			#store xmm3 to forth 16 bytes from second argument


	movq   0(%rsp), %r11
	movq   %r11, %rsp

ret

