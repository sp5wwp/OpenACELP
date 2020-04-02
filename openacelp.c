//Standard includes
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#define DEBUG										//comment it out later

//-----------------------------ACELP includes & defines------------------------------
//Global consts
#include "N_POW2.h"									//powers of 2 array
#include "LOG2.h"									//array for Log2 computation

#define FRAME_LEN		20							//voice frame length in ms
#define FRAME_SIZ		160							//voice frame samples number, 0.02s*8000Hz
#define	WINDOW_LEN		30							//length of window for autocorrelation computation
#define WINDOW_SIZ		240							//window samples number, 0.03s*8000Hz
#define ALPHA			(float)32735.0/32768.0		//alpha coeff for the pre-processing filter

//Global vars
float    w[WINDOW_SIZ];								//modified Hamming window w(n) coeffs for speech analysis
uint32_t r[11];										//autocorrelation values

//Global flags
//

//----------------------------------ACELP functions----------------------------------
//input speech signal pre-processing
//y[i] = x[i]/2 - x[i-1]/2 + alpha * y[i-1]
//arg1: present frame, arg2: output
void Speech_Pre_Process(int16_t *inp, int16_t *outp)
{
	#ifdef DEBUG
	if(inp==NULL || outp==NULL)
	{
		printf("NULL pointer at Speech_Pre_Process()\n");
		exit(0);
	}
	#endif
	
	outp[0]=inp[0];
	
	for(uint8_t i=1; i<WINDOW_SIZ; i++)
	{
		outp[i] = inp[i]/2 - inp[i-1]/2 + ALPHA*outp[i-1];
	}
}

//compute the modified Hamming window w(n) coeffs for speech analysis
void Analysis_Window_Init(float *w)
{
	uint8_t L2 = 40;			//40 samples look ahead
	uint8_t L1 = WINDOW_SIZ-40;
	
	for(uint8_t i=0; i<L1; i++)
	{
		w[i] = 0.54 - 0.46 * cos((M_PI*i)/((float)L1-1.0));
	}
	for(uint8_t i=L1; i< L1+L2; i++)
	{
		w[i] = 0.54 + 0.46 * cos((M_PI*(i-L1))/((float)L2-1.0));
	}
}

//multiply processed speech samples with modified Hamming window
void Window_Speech(int16_t *inp, int16_t *outp)
{
	for(uint8_t i=0; i<WINDOW_SIZ; i++)
		outp[i] = inp[i] * w[i];
}

//autocorrelation r(k) computation, k=0..10
void Autocorr(int16_t *spch, uint32_t *r)
{
	uint8_t ovf;
	uint32_t prev;
	uint8_t norm_shift=0;		//shifts left needed to normalize r[]

	memset((uint8_t*)r, 0, 11*sizeof(uint32_t));
	
	//initially, set r[0]=1 to avoid r[] containing only zeros
	r[0]=1;
	
	//r[0] calculation
	//if r[0] overflows uint32_t, divide the signal by 4
	do
	{
		ovf = 0;
		
		for(uint8_t i=0; i<WINDOW_SIZ; i++)
		{
			prev = r[0];
			r[0] += (uint32_t)(spch[i] * spch[i]);
		
			if(r[0] < prev)	//overflow occured?
			{				
				//divide the signal by 4
				for(uint8_t j=0; j<WINDOW_SIZ; j++)
					spch[j]/=4;
				
				//zero r[]
				memset((uint8_t*)r, 0, 11*sizeof(uint32_t));
				
				ovf = 1;
				
				//break the "for" loop
				break;
			}
		}
	}
	while(ovf);
	
	#ifdef DEBUG
	printf("r[0]=%lld\n", r[0]);
	#endif

	//r[0] normalization to the uint32_t limit
	//shift left until the MSB==1
	for(uint8_t i=0; i<32; i++)
	{
		while(!(r[0] & 0x80000000))
		{
			r[0]<<=1;
			norm_shift++;
		}
	}
	
	#ifdef DEBUG
	printf("norm_shift=%d\n", norm_shift);
	#endif
}

//LPC coefficients calculation
//based on the Levison-Durbin algorithm
void LD_Solver(uint32_t *r, int16_t *A)
{
	;
}

//innovative codebook search
int16_t Innov_Codebook_Search(int16_t dn[], int16_t f[], int16_t h[], int16_t rr[][32],
             int16_t cod[], int16_t y[], int16_t *sign, int16_t *shift_code)
{
	;
}

int main(void)
{
	int16_t tst_spch[WINDOW_SIZ];
	int16_t tst_out[WINDOW_SIZ];
	
	for(uint8_t i=0; i<WINDOW_SIZ; i++)
		tst_spch[i]=(sin(i/8.0)*0x1FFF)*2;
	
	Speech_Pre_Process(tst_spch, tst_out);
	
	for(uint8_t i=0; i<WINDOW_SIZ; i++)
		;//printf("%d|%d\n", tst_spch[i], tst_out[i]);
	
	Analysis_Window_Init(w);
	
	for(uint8_t i=0; i<WINDOW_SIZ; i++)
		;//printf("%f\n", w[i]);
		
	memcpy((int16_t*)tst_spch, (int16_t*)tst_out, WINDOW_SIZ*sizeof(int16_t));
		
	Window_Speech(tst_spch, tst_out);
	
	for(uint8_t i=0; i<WINDOW_SIZ; i++)
		;//printf("%d\n", tst_out[i]);
		
	Autocorr(tst_out, r);
	
	return 0;
}

