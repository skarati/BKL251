.p2align 5
.globl conditionalSwap
conditionalSwap:

	movdqa	0(%rdi),%xmm1			#store first 16 bytes from second argument to xmm1
	movdqa	16(%rdi),%xmm2			#store second 16 bytes from second argument to xmm2
	movdqa	32(%rdi),%xmm3			#store third 16 bytes from second argument to xmm3
	movdqa	48(%rdi),%xmm4			#store forth 16 bytes from second argument to xmm4

	movdqa	0(%rsi),%xmm5			#store first 16 bytes from third argument to xmm5
	movdqa	16(%rsi),%xmm6			#store second 16 bytes from third argument to xmm6
	movdqa	32(%rsi),%xmm7			#store third 16 bytes from third argument to xmm7
	movdqa	48(%rsi),%xmm8			#store forth 16 bytes from third argument to xmm8

	movdqa %xmm5, %xmm9
	movdqa %xmm6, %xmm10
	movdqa %xmm7, %xmm11
	movdqa %xmm8, %xmm12

	pxor %xmm1, %xmm9
	pxor %xmm2, %xmm10
	pxor %xmm3, %xmm11
	pxor %xmm4, %xmm12

	pand %xmm0, %xmm9
	pand %xmm0, %xmm10
	pand %xmm0, %xmm11
	pand %xmm0, %xmm12

	pxor %xmm9, %xmm1
	pxor %xmm10, %xmm2
	pxor %xmm11, %xmm3
	pxor %xmm12, %xmm4

	pxor %xmm9, %xmm5
	pxor %xmm10, %xmm6
	pxor %xmm11, %xmm7
	pxor %xmm12, %xmm8

	movdqa	%xmm1,0(%rdi)
	movdqa	%xmm2,16(%rdi)
	movdqa	%xmm3,32(%rdi)
	movdqa	%xmm4,48(%rdi)

	movdqa	%xmm5,0(%rsi)
	movdqa	%xmm6,16(%rsi)
	movdqa	%xmm7,32(%rsi)
	movdqa	%xmm8,48(%rsi)

	movdqa	0(%rdx),%xmm1			#store first 16 bytes from second argument to xmm1
	movdqa	16(%rdx),%xmm2			#store second 16 bytes from second argument to xmm2
	movdqa	32(%rdx),%xmm3			#store third 16 bytes from second argument to xmm3
	movdqa	48(%rdx),%xmm4			#store forth 16 bytes from second argument to xmm4

	movdqa	0(%rcx),%xmm5			#store first 16 bytes from third argument to xmm5
	movdqa	16(%rcx),%xmm6			#store second 16 bytes from third argument to xmm6
	movdqa	32(%rcx),%xmm7			#store third 16 bytes from third argument to xmm7
	movdqa	48(%rcx),%xmm8			#store forth 16 bytes from third argument to xmm8

	movdqa %xmm5, %xmm9
	movdqa %xmm6, %xmm10
	movdqa %xmm7, %xmm11
	movdqa %xmm8, %xmm12

	pxor %xmm1, %xmm9
	pxor %xmm2, %xmm10
	pxor %xmm3, %xmm11
	pxor %xmm4, %xmm12

	pand %xmm0, %xmm9
	pand %xmm0, %xmm10
	pand %xmm0, %xmm11
	pand %xmm0, %xmm12

	pxor %xmm9, %xmm1
	pxor %xmm10, %xmm2
	pxor %xmm11, %xmm3
	pxor %xmm12, %xmm4

	pxor %xmm9, %xmm5
	pxor %xmm10, %xmm6
	pxor %xmm11, %xmm7
	pxor %xmm12, %xmm8

	movdqa	%xmm1,0(%rdx)
	movdqa	%xmm2,16(%rdx)
	movdqa	%xmm3,32(%rdx)
	movdqa	%xmm4,48(%rdx)

	movdqa	%xmm5,0(%rcx)
	movdqa	%xmm6,16(%rcx)
	movdqa	%xmm7,32(%rcx)
	movdqa	%xmm8,48(%rcx)

ret
