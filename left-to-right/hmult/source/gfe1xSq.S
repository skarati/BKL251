.p2align 5
.globl gfe1xnSq
gfe1xnSq:

	movdqa	0(%rsi),%xmm0			#store first 16 bytes from second argument to xmm0
	movdqa	16(%rsi),%xmm2			#store second 16 bytes from second argument to xmm2
	movdqa	32(%rsi),%xmm4			#store third 16 bytes from second argument to xmm4
	movdqa	48(%rsi),%xmm6			#store forth 16 bytes from second argument to xmm6

	movq    %rdx, %r8			#r8 = rdx (third argument that is the counter)

    	.LOOP:					#start the loop

		subq    $1, %r8				# reduce the counter by 1

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

		cmpq    $0, %r8			#compare whether the counter is 0 or not

    	jne .LOOP				#if counter is not zero, jump to label LOOP


	movdqa	%xmm0,0(%rdi)			#store xmm0 to first 16 bytes from first argument
	movdqa	%xmm2,16(%rdi)			#store xmm2 to second 16 bytes from first argument
	movdqa	%xmm4,32(%rdi)			#store xmm4 to third 16 bytes from first argument
	movdqa	%xmm6,48(%rdi)			#store xmm6 to forth 16 bytes from first argument

ret
