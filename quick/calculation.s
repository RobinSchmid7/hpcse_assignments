	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 10, 15	sdk_version 10, 15
	.section	__TEXT,__literal8,8byte_literals
	.p2align	3               ## -- Begin function main
LCPI0_0:
	.quad	4607182418800017408     ## double 1
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_main
	.p2align	4, 0x90
_main:                                  ## @main
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -24
	movq	%rsi, %rbx
	cmpl	$2, %edi
	jg	LBB0_2
## %bb.1:
	movq	(%rbx), %rsi
	leaq	L_.str(%rip), %rdi
	xorl	%eax, %eax
	callq	_printf
	jmp	LBB0_3
LBB0_2:
	movq	8(%rbx), %rdi
	callq	_atof
	movsd	%xmm0, -16(%rbp)        ## 8-byte Spill
	movq	16(%rbx), %rdi
	callq	_atof
	movsd	LCPI0_0(%rip), %xmm1    ## xmm1 = mem[0],zero
	movapd	%xmm1, %xmm2
	movsd	-16(%rbp), %xmm3        ## 8-byte Reload
                                        ## xmm3 = mem[0],zero
	subsd	%xmm3, %xmm2
	divsd	%xmm0, %xmm3
	addsd	%xmm2, %xmm3
	divsd	%xmm3, %xmm1
	leaq	L_.str.1(%rip), %rdi
	movapd	%xmm1, %xmm0
	movb	$1, %al
	callq	_printf
LBB0_3:
	xorl	%eax, %eax
	addq	$8, %rsp
	popq	%rbx
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__cstring,cstring_literals
L_.str:                                 ## @.str
	.asciz	"Usage: %s p n \n"

L_.str.1:                               ## @.str.1
	.asciz	"s = %f\n"


.subsections_via_symbols
