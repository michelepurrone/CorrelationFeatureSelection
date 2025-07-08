;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"
        
section .data ; Sezione contenente dati inizializzati

    align 16
    uno: dd 1.0, 1.0, 1.0, 1.0
    
section .bss ; Sezione contenente dati non inizializzati
    
    alignb 16
        temp resd 0

    alignb 16
        i	resd 	1
    alignb 16
        j	resd	1
    alignb 16
        nx	resd	1

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------
        
; ------------------------------------------------------------ CORRELATION_CF ------------------------------------------------------------
;extern void correlation_cf(type* feature, type* labels, int n, type* r_cf);
section .text

global correlation_cf
	;feature		equ	8
	;labels		    equ	12
	;n		        equ	16
	;r_cf		    equ	20
	
correlation_cf:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp
		mov		ebp, esp
		push		ebx
		push		esi
		push		edi

		; ------------------------------------------------------------
		; FUNZIONE ASSEMBLY
		; ------------------------------------------------------------
		mov	eax, [ebp+8]		;feature addr
		mov	ebx, [ebp+12]		;labels addr
		mov	ecx, [ebp+16]		;n
		
		xor	esi, esi		;index

		xorps	xmm0, xmm0		;sum 0
		xorps	xmm1, xmm1		;sum 1
		xorps	xmm2, xmm2		;mu 0
		xorps	xmm3, xmm3		;mu 1
		xorps	xmm4, xmm4		

		xorps xmm7, xmm7 ;vettore 0
		
		xor edx, edx ;n_0
		xor edi, edi ;n_1

ciclo_if01:
		add	esi, 1
		cmp	esi, ecx ;for(...i < n;...)
		jg	fine_ciclo_if01
		
		movss	xmm5, [eax+4*esi-4]	;feature[i]
		movss	xmm6, [ebx+4*esi-4]	;labels[i]
		
		ucomiss xmm6, xmm7 ;if(labels[i] == 0)
		je increment_n0
		
		jmp increment_n1
		
increment_n0:
        addss xmm0, xmm5 ;sum_0 += feature[i]
        inc edx ;n_0++
        jmp	ciclo_if01

increment_n1: ;else		
		addss xmm1, xmm5 ;sum_1 += feature[i]
		inc edi ;n_1++
		jmp	ciclo_if01


fine_ciclo_if01:
        movss xmm2, xmm0
        cvtsi2ss xmm4, edx		;n_0 replicato 4 volte
		divss xmm2, xmm4 ;mu_0 = sum_0 / n_0
		
		movss xmm3, xmm1
		cvtsi2ss xmm5, edi		;n_1 replicato 4 volte
		divss xmm3, xmm5 ;mu_1 = sum_1 / n_1
		
		;type mu = (sum_0 + sum_1) / n
		cvtsi2ss xmm7, ecx		;n replicato 4 volte
		addss xmm0, xmm1 ;sum_0 + sum_1
		divss xmm0, xmm7 ;mu = (sum_0 + sum_1) / n
		
		movaps xmm1, [uno]
		subss xmm7, xmm1 ;(n - 1.0)
		divss xmm1, xmm7 ;inv_n_minus_1 = 1.0 / (n - 1.0)
		
		cvtsi2ss xmm7, ecx		;n replicato 4 volte
		mulss xmm7, xmm7 ;n * n
		mulss xmm4, xmm5 ;n_0 * n_1
		divss xmm4, xmm7 ;(n_0 * n_1) / (n * n)
		sqrtss xmm4, xmm4 ;sqrt_n0_n1_n_n = sqrt((n_0 * n_1) / (type) (n * n))
		
		shufps xmm0, xmm0, 00000000b
		xor	esi, esi		;index
		xorps xmm6, xmm6

ciclo_sommatoria_SIMD:
		add	esi, 4
		cmp	esi, ecx ;for(...i < n;...)
		jg	fine_ciclo_sommatoria_SIMD
		
		movaps	xmm5, [eax+4*esi-16]	;feature[i]
		
		
		subps xmm5, xmm0 ;(feature[i] - mu)
		
		mulps xmm5, xmm5 ;(feature[i] - mu) * (feature[i] - mu)
		
		addps xmm6, xmm5 ;sommatoria += (feature[i] - mu) * (feature[i] - mu)
		
		jmp ciclo_sommatoria_SIMD
		
fine_ciclo_sommatoria_SIMD:
        sub esi, 4
        
ciclo_sommatoria:
        add	esi, 1
		cmp	esi, ecx ;for(...i < n;...)
		jg	fine_ciclo_sommatoria
		
		movss	xmm5, [eax+4*esi-4]	;feature[i]
		subss xmm5, xmm0 ;(feature[i] - mu)
		mulss xmm5, xmm5 ;(feature[i] - mu) * (feature[i] - mu)
		
		addss xmm6, xmm5 ;sommatoria += (feature[i] - mu) * (feature[i] - mu)
		
		jmp ciclo_sommatoria
		
fine_ciclo_sommatoria:

        haddps	xmm6, xmm6
		haddps	xmm6, xmm6		;prima cella sommatoria
		
        mulss xmm1, xmm6 ;inv_n_minus_1 * sommatoria
        sqrtss xmm1, xmm1 ;dev = sqrt(inv_n_minus_1 * sommatoria)
        
        subss xmm2, xmm3 ;(mu_0 - mu_1)
        divss xmm2, xmm1 ;((mu_0 - mu_1) / dev)
        mulss xmm2, xmm4 ;r_cf = ((mu_0 - mu_1) / dev) * sqrt_n0_n1_n_n
		
		mov	edx, [ebp+20]
		movss	[edx], xmm2		;ret corr value

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret	

; ------------------------------------------------------------ CORRELATION_CF ------------------------------------------------------------

; ------------------------------------------------------------ CORRELATION_FF ------------------------------------------------------------
;extern void correlation_ff(type* x, type* y, int n, type* r_ff);
section .text

global correlation_ff
	;x		    equ	8
	;y		    equ	12
	;n		    equ	16
	;r_ff		equ	20
	
correlation_ff:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp
		mov		ebp, esp
		push		ebx
		push		esi
		push		edi

		; ------------------------------------------------------------
		; FUNZIONE ASSEMBLY
		; ------------------------------------------------------------
		mov	eax, [ebp+8]		;x addr
		mov	ebx, [ebp+12]		;y addr
		mov	ecx, [ebp+16]		;n
		
		xor	esi, esi		;index

		xorps	xmm0, xmm0		;sum X
		xorps	xmm1, xmm1		;sum Y
		xorps	xmm2, xmm2		;sum X*Y
		xorps	xmm3, xmm3		;square X
		xorps	xmm4, xmm4		;square Y

ciclo_SIMD:
		add	esi, 4
		cmp	esi, ecx
		jg	fine_ciclo_SIMD
		movaps	xmm5, [eax+4*esi-16]	;x[i]
		movaps	xmm6, [ebx+4*esi-16]	;y[i]
		addps	xmm0, xmm5		;sumX
		addps	xmm1, xmm6		;sumY
		movaps	xmm7, xmm5		;copy x[i]
		mulps	xmm7, xmm6		;x[i]*y[i]
		addps	xmm2, xmm7		;sumX*Y
		mulps	xmm5, xmm5		;x[i]^2
		mulps	xmm6, xmm6		;y[i]^2
		addps	xmm3, xmm5		;sumX^2
		addps	xmm4, xmm6		;sumY^2

		jmp	ciclo_SIMD

fine_ciclo_SIMD:
		sub	esi, 4

ciclo_1:
		add	esi, 1
		cmp	esi, ecx
		jg	fine_ciclo
		movss	xmm5, [eax+4*esi-4]	;x[i]
		movss	xmm6, [ebx+4*esi-4]	;y[i]
		addss	xmm0, xmm5		;sumX
		addss	xmm1, xmm6		;sumY
		movss	xmm7, xmm5		;copy x[i]
		mulss	xmm7, xmm6		;x[i]*y[i]
		addss	xmm2, xmm7		;sumX*Y
		mulss	xmm5, xmm5		;x[i]^2
		mulss	xmm6, xmm6		;y[i]^2
		addss	xmm3, xmm5		;sumX^2
		addss	xmm4, xmm6		;sumY^2

		jmp	ciclo_1

fine_ciclo:
		haddps	xmm0, xmm0
		haddps	xmm0, xmm0		;prima cella sumX
		haddps	xmm1, xmm1
		haddps	xmm1, xmm1		;prima cella sumY
		haddps	xmm2, xmm2
		haddps	xmm2, xmm2		;prima cella sumX*Y
		haddps	xmm3, xmm3
		haddps	xmm3, xmm3		;prima cella sumX^2
		haddps	xmm4, xmm4
		haddps	xmm4, xmm4		;prima cella sumY^2

		xorps	xmm5, xmm5
		xorps	xmm6, xmm6
		xorps	xmm7, xmm7

		;( N*xmm2 - xmm0 * xmm1 ) / sqrt[ ( N*xmm3 - xmm0*xmm0) * ( N*xmm4 - xmm1*xmm1 ) ]
		;    V1         V2			V3	V4		V5	V6
		;'_______________________'      '______________________'  '______________________'
		;	P1				P2				P3

		cvtsi2ss xmm5, ecx		;n replicato 4 volte
		mulss	xmm2, xmm5		;V1
		movss	xmm6, xmm0
		mulss	xmm6, xmm1		;V2
		mulss	xmm3, xmm5		;V3
		mulss	xmm0, xmm0		;V4
		mulss	xmm4, xmm5		;V5
		mulss	xmm1, xmm1		;V6

		subss	xmm2, xmm6		;P1
		subss	xmm3, xmm0		;P2
		subss	xmm4, xmm1		;P3

		mulss	xmm3, xmm4		;contenuto sqrt
		sqrtss	xmm3, xmm3		;sqrt
		
		divss	xmm2, xmm3		;corr value
		
		mov	edx, [ebp+20]
		movss	[edx], xmm2		;ret corr value

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret	
; ------------------------------------------------------------ CORRELATION_FF ------------------------------------------------------------
    
; ------------------------------------------------------------ TRASPOSTA ------------------------------------------------------------
;extern void transposeOMP(MATRIX A, MATRIX B, int ncols, int nrows, int n);
section .text

        A	        equ	8
        B	        equ 12
        c	        equ	16
        r	        equ	20
        n           equ 24

global transposeOMP

transposeOMP:
        ; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp
		mov		ebp, esp
		push		ebx
		push		esi
		push		edi
		
		; ------------------------------------------------------------
		; FUNZIONE ASSEMBLY
		; ------------------------------------------------------------
        mov     eax, [ebp+n]
        mov     ecx, [ebp+r]
        cdq
        div     ecx                 ; n/indiceRiga
        mov     esi, eax            ; i
        mov     edi, edx            ; j
        mov     eax, [ebp+c]
        mul     edi                 ; eax = indiceCol*j
        add     eax, esi            ; eax = indiceCol*j+i
        mov     ebx, [ebp+A]
        movss   xmm0, [ebx+eax*4]       ; xmm0 = k[indiceCol*j+i]
        mov     ecx, [ebp+B]
        mov     edx, [ebp+n]
        movss   [ecx+edx*4], xmm0       ; k[n] = xmm0

        ; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret	

; ------------------------------------------------------------ TRASPOSTA ------------------------------------------------------------
