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
;     nasm -f elf64 regression64.nasm
;
%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
    align 32
    uno dq 1.0, 1.0, 1.0, 1.0


section .bss			; Sezione contenente dati non inizializzati
    alignb 32
    s resq 4

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

; ------------------------------------------------------------ CORRELATION_CF ------------------------------------------------------------
;extern void correlation_cf(type* feature, type* labels, int n, type* r_cf);
section .text
 
global correlation_cf
	;feature		rdi
	;labels			rsi
	;n			    rdx
	;r_cf			rcx
	
correlation_cf:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp
		mov		rbp, rsp
		pushaq
 
		; ------------------------------------------------------------
		; FUNZIONE ASSEMBLY
		; ------------------------------------------------------------
 
		vxorpd	ymm0, ymm0, ymm0		;sum_0
		vxorpd	ymm1, ymm1, ymm1		;sum_1
		vxorpd	ymm2, ymm2, ymm2		;mu_0
		vxorpd	ymm3, ymm3, ymm3		;mu_1
		vxorpd	ymm4, ymm4, ymm4		;vettore 0
		
		xor r10, r10 ;n_0
		xor r11, r11 ;n_1
		
		xor	rax, rax		;index
 
ciclo_if01:
		add	rax, 1
		cmp	rax, rdx ;for(...i < n;...)
		jg	fine_ciclo_if01
		
		vmovsd	xmm5, [rdi+8*rax-8]	;feature[i]
		vmovsd	xmm6, [rsi+8*rax-8]	;labels[i]
		
		vcomisd xmm6, xmm4
		je	increment_n0
		jmp	increment_n1
		
increment_n0:
	    vaddsd	xmm0, xmm5	;sum_0 += feature[i]
	    add	r10, 1		;n_0++
	    jmp	ciclo_if01
 
increment_n1: ;else		
		vaddsd	xmm1, xmm5 	;sum_1 += feature[i]
		add	r11, 1		;n_1++
		jmp	ciclo_if01
 
 
fine_ciclo_if01:
        
		vcvtsi2sd xmm5, xmm5, r10		;n_0 replicato
		vcvtsi2sd xmm6, xmm6, r11		;n_1 replicato
		vcvtsi2sd xmm7, xmm7, rdx		;n replicato
 
		vmovsd xmm2, xmm0
		vdivsd xmm2, xmm5 	;mu_0 = sum_0 / n_0
		
		vmovsd xmm3, xmm1
		vdivsd xmm3, xmm6 	;mu_1 = sum_1 / n_1
		
		;type mu = (sum_0 + sum_1) / n;
		vaddsd xmm0, xmm1 	;sum_0 + sum_1
		vdivsd xmm0, xmm7 	;mu = (sum_0 + sum_1) / n
		
		vmovsd xmm1, [uno] 
		vmovsd xmm4, xmm7	;n
		vsubsd xmm4, xmm1 	;(n - 1.0)
		vdivsd xmm1, xmm4 	;inv_n_minus_1 = 1.0 / (n - 1.0)
		
		vmulsd xmm7, xmm7 	;n * n
		vmulsd xmm5, xmm6 	;n_0 * n_1
		vdivsd xmm5, xmm7 	;(n_0 * n_1) / (n * n)
		vsqrtsd xmm5, xmm5 	;sqrt_n0_n1_n_n = sqrt((n_0 * n_1) / (type) (n * n))
		
		 
		vbroadcastsd	ymm0, xmm0
		vxorpd	ymm4, ymm4, ymm4
		vxorpd	ymm6, ymm6, ymm6	;sommatoria
		xor	rax, rax		;index
 
ciclo_sommatoriaSIMD:
		add	rax, 4
		cmp	rax, rdx ;for(...i < n;...)
		jg	fine_ciclo_sommatoriaSIMD
		
		vmovapd	ymm4, [rdi+8*rax-32]	;feature[i]
		vsubpd	ymm4, ymm0 		;(feature[i] - mu)
		vmulpd	ymm4, ymm4 		;(feature[i] - mu) * (feature[i] - mu)
		vaddpd 	ymm6, ymm4 		;sommatoria += (feature[i] - mu) * (feature[i] - mu);
			
		jmp ciclo_sommatoriaSIMD
		
fine_ciclo_sommatoriaSIMD:
		sub	rax, 4
		
		vxorpd ymm8, ymm8, ymm8
        vxorpd ymm9, ymm9, ymm9  
        vhaddpd ymm8, ymm6, ymm6  
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm6, xmm9, xmm8
	
 
ciclo_sommatoria:
		add	rax, 1
		cmp	rax, rdx ;for(...i < n;...)
		jg	fine_ciclo_sommatoria
		
		vmovsd	xmm4, [rdi+8*rax-8]	;feature[i]
		vsubsd	xmm4, xmm0 		;(feature[i] - mu)
		vmulsd	xmm4, xmm4 		;(feature[i] - mu) * (feature[i] - mu)	
		vaddsd 	xmm6, xmm4 		;sommatoria += (feature[i] - mu) * (feature[i] - mu)
		
		jmp ciclo_sommatoria
		
fine_ciclo_sommatoria:
        vmulsd	xmm1, xmm6 		;inv_n_minus_1 * sommatoria
        vsqrtsd	xmm1, xmm1 		;dev = sqrt(inv_n_minus_1 * sommatoria)
        
        vsubsd	xmm2, xmm3 		;(mu_0 - mu_1)
        vdivsd	xmm2, xmm1 		;((mu_0 - mu_1) / dev)
       	vmulsd	xmm2, xmm5 		;r_cf = ((mu_0 - mu_1) / dev) * sqrt_n0_n1_n_n
       	
		vmovsd   [rcx], xmm2         ;ret corr value
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq
		mov		rsp, rbp
		pop		rbp
		ret
; ------------------------------------------------------------ CORRELATION_CF ------------------------------------------------------------
 
; ------------------------------------------------------------ CORRELATION_FF ------------------------------------------------------------
;extern void correlation_ff(type* x, type* y, int n, type* r_ff);
section .text

global correlation_ff

correlation_ff:
        ; ------------------------------------------------------------
	    ; Sequenza di ingresso nella funzione
	    ; ------------------------------------------------------------
		push	rbp
		mov		rbp, rsp
		pushaq
		
		; ------------------------------------------------------------
		; FUNZIONE ASSEMBLY
		; ------------------------------------------------------------
		
        vxorpd  ymm0, ymm0, ymm0  ; sum X
        vxorpd  ymm1, ymm1, ymm1  ; sum Y
        vxorpd  ymm2, ymm2, ymm2  ; sum X*Y
        vxorpd  ymm3, ymm3, ymm3  ; square X
        vxorpd  ymm4, ymm4, ymm4  ; square Y
    
        xor     rax, rax    ; index

ciclo_SIMD:

        add     rax, 4
        cmp     rax, rdx
        jg      fine_ciclo_SIMD
        vmovapd ymm5, [rdi+8*rax-32]  ; x[i]
        vmovapd ymm6, [rsi+8*rax-32]  ; y[i]
        vaddpd  ymm0, ymm0, ymm5      ; sumX
        vaddpd  ymm1, ymm1, ymm6      ; sumY
        vmovapd ymm7, ymm5            ; copy x[i]
        vmulpd  ymm7, ymm7, ymm6      ; x[i]*y[i]
        vaddpd  ymm2, ymm2, ymm7      ; sumX*Y
        vmulpd  ymm5, ymm5, ymm5      ; x[i]^2
        vmulpd  ymm6, ymm6, ymm6      ; y[i]^2
        vaddpd  ymm3, ymm3, ymm5      ; sumX^2
        vaddpd  ymm4, ymm4, ymm6      ; sumY^2

        jmp     ciclo_SIMD

fine_ciclo_SIMD:

        sub     rax, 4
        
        vxorpd ymm8, ymm8, ymm8
        vxorpd ymm9, ymm9, ymm9  
        vhaddpd ymm8, ymm0, ymm0  
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm0, xmm9, xmm8
 
        vhaddpd ymm8, ymm1, ymm1   
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm1, xmm9, xmm8
         
        vhaddpd ymm8, ymm2, ymm2   
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm2, xmm9, xmm8
         
        vhaddpd ymm8, ymm3, ymm3   
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm3, xmm9, xmm8
         
        vhaddpd ymm8, ymm4, ymm4   
        vextractf128 xmm9, ymm8, 1
        vaddsd xmm4, xmm9, xmm8               


ciclo_1:
        add     rax, 1
        cmp     rax, rdx
        jg      fine_ciclo
        vmovsd  xmm5, [rdi+8*rax-8]  ; x[i]
        vmovsd  xmm6, [rsi+8*rax-8]  ; y[i]
        vaddsd  xmm0, xmm0, xmm5     ; sumX
        vaddsd  xmm1, xmm1, xmm6     ; sumY
        vmovsd  xmm7, xmm5           ; copy x[i]
        vmulsd  xmm7, xmm7, xmm6     ; x[i]*y[i]
        vaddsd  xmm2, xmm2, xmm7     ; sumX*Y
        vmulsd  xmm5, xmm5, xmm5     ; x[i]^2
        vmulsd  xmm6, xmm6, xmm6     ; y[i]^2
        vaddsd  xmm3, xmm3, xmm5     ; sumX^2
        vaddsd  xmm4, xmm4, xmm6     ; sumY^2

        jmp     ciclo_1

fine_ciclo:
		
        vxorpd  ymm5, ymm5, ymm5
        vxorpd  ymm6, ymm6, ymm6
        vxorpd  ymm7, ymm7, ymm7
        
        ;( N*xmm2 - xmm0 * xmm1 ) / sqrt[ ( N*xmm3 - xmm0*xmm0) * ( N*xmm4 - xmm1*xmm1 ) ]
		;    V1         V2			V3	V4		V5	V6
		;'_______________________'      '______________________'  '______________________'
		;	P1				P2				P3

        ; Calcolo della correlazione
        vcvtsi2sd xmm5, xmm5, rdx    ; n replicato 4 volte
        vmulsd   xmm2, xmm2, xmm5    ; V1
        vmovsd   xmm6, xmm0
        vmulsd   xmm6, xmm6, xmm1    ; V2
        vmulsd   xmm3, xmm3, xmm5    ; V3
        vmulsd   xmm0, xmm0, xmm0    ; V4
        vmulsd   xmm4, xmm4, xmm5    ; V5
        vmulsd   xmm1, xmm1, xmm1    ; V6

        vsubsd   xmm2, xmm2, xmm6    ; P1
        vsubsd   xmm3, xmm3, xmm0    ; P2
        vsubsd   xmm4, xmm4, xmm1    ; P3

        vmulsd   xmm3, xmm3, xmm4    ; contenuto sqrt
        vsqrtsd  xmm3, xmm3, xmm3    ; sqrt
    
        vdivsd   xmm2, xmm2, xmm3    ; corr value

        vmovsd   [rcx], xmm2         ; ret corr value

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq
		mov		rsp, rbp
		pop		rbp
		ret
		
; ------------------------------------------------------------ CORRELATION_FF ------------------------------------------------------------
		
; ------------------------------------------------------------ TRASPOSTA ------------------------------------------------------------

section .text
        ;A       rdi
        ;B       rsi
        ;col     rdx
        ;rig     rcx

global transpose

transpose:

        ; ------------------------------------------------------------
	    ; Sequenza di ingresso nella funzione
	    ; ------------------------------------------------------------
		push	rbp
		mov		rbp, rsp
		pushaq
        
        ; ------------------------------------------------------------
        ; FUNZIONE ASSEMBLY
        ; ------------------------------------------------------------

        xor         r12, r12            ; n = 0
        mov         r13, rdx            ; r13 = col
        mov         rax, r13            ; rax = col
        mul         rcx                 ; rax *= rig
        mov         r9, rax             ; r9 = col*rig
ciclot: cmp         r12, r9             ; n < col*rig ?
        jge         finet               ; se <= salta a finet
        mov         rax, r12            ; rax = n
        div         rcx                 ; rax = n / rig, rdx = n % rig
        mov         r10, rax            ; r10 = i
        mov         r11, rdx            ; r11 = j
        mov         rax, r13            ; rax = col
        mul         r11                 ; rax = col*j
        add         rax, r10            ; rax += i
        vmovsd      xmm1, [rdi+rax*8]   
        vmovsd      [rsi+r12*8], xmm1   
        inc         r12                 ; n++
        jmp         ciclot              ; salta ciclot

finet:  
        ; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq
		mov		rsp, rbp
		pop		rbp
		ret
                
; ------------------------------------------------------------ TRASPOSTA ------------------------------------------------------------

