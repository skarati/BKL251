.p2align 5
.globl gfe1xAdd
gfe1xAdd:

	movdqa	0(%rsi),%xmm0			#store first 16 bytes from second argument to xmm0
	movdqa	16(%rsi),%xmm1			#store second 16 bytes from second argument to xmm1
	movdqa	32(%rsi),%xmm2			#store third 16 bytes from second argument to xmm2
	movdqa	48(%rsi),%xmm3			#store forth 16 bytes from second argument to xmm3

	movdqa	0(%rdx),%xmm4			#store first 16 bytes from third argument to xmm4
	movdqa	16(%rdx),%xmm5			#store second 16 bytes from third argument to xmm5
	movdqa	32(%rdx),%xmm6			#store third 16 bytes from third argument to xmm6
	movdqa	48(%rdx),%xmm7			#store forth 16 bytes from third argument to xmm7

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4 = xmm0 + xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5 = xmm1 + xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6 = xmm2 + xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa	%xmm0,0(%rdi)			#store xmm0 to first 16 bytes from first argument
	movdqa	%xmm1,16(%rdi)			#store xmm1 to second 16 bytes from first argument
	movdqa	%xmm2,32(%rdi)			#store xmm2 to third 16 bytes from first argument
	movdqa	%xmm3,48(%rdi)			#store xmm3 to forth 16 bytes from first argument

ret
