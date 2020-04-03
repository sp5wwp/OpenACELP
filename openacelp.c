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
//additional bandwidth expansion f=60Hz for f_s=8000Hz sample rate
void Autocorr(int16_t *spch, uint32_t *r)
{
	uint8_t ovf;
	uint32_t prev;
	uint8_t norm_shift=0;		//shifts left needed to normalize r[]

	//initially, set r[0]=1 to avoid r[] containing only zeros
	//zero out the rest
	r[0]=1;
	memset((uint8_t*)&r[1], 0, 10*sizeof(uint32_t));
	
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
					spch[j] /= 4;
				
				//re-set r[0] to 1
				r[0] = 1;
				
				ovf = 1;
				
				//break the "for" loop
				break;
			}
		}
	}
	while(ovf);

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
	
	//r[1]..r[10] calculation
	for(uint8_t i=1; i<=10; i++)
	{
		for(uint8_t j=0; j<WINDOW_SIZ; j++)
			r[i] += spch[j] * spch[j-i];
			
		//normalize
		r[i] <<= norm_shift;
	}
	
	#ifdef DEBUG
	printf("norm_shift=%d\n", norm_shift);
	for(uint8_t i=0; i<=10; i++)
		printf("r[%d]=%lld\n", i, r[i]);
	printf("\n");
	#endif
	
	//bandwidth expansion, f=60Hz, f_s=8000Hz
	//r[0] is multiplied by 1.00005, which is equivalent to adding a noise floor at -43 dB
	float w_lag[11]={1.0};	//window for bandwidth expansion, w_lag[0]=1.0
	
	for(uint8_t i=1; i<=10; i++)
	{
		w_lag[i] = exp(-0.5 * (2.0 * M_PI * 60.0 * i)/(8000.0));
		w_lag[i] /= 1.00005;
	}
	
	#ifdef DEBUG
	for(uint8_t i=0; i<=10; i++)
		printf("w_lag[i]=%f\n", w_lag[i]);
	printf("\n");
	#endif
	
	//window the autocorrelation values
	for(uint8_t i=0; i<11; i++)
		r[i] *= w_lag[i];
	
	#ifdef DEBUG
	for(uint8_t i=0; i<=10; i++)
		printf("r[%d]=%lld\n", i, r[i]);
	printf("\n");
	#endif
}

//LP coefficients calculation
//based on the Levison-Durbin algorithm
//arg1: modified autocorrelation matrix, arg2: LP filter coeffs
void LD_Solver(uint32_t *r, int32_t *a)
{
	float v, c, e;
	float rv[11];	//
	float tmp[11];	//
	//float E[11];
	//float sum;
	
	/*E[0] = (float)r[0];
	
	for(uint8_t i=1; i<=10; i++)
	{
		sum=0.0;
		
		if(i>1)
		{
			sum=0.0;
			for(uint8_t j=1; j<=i-1; j++)
				sum += (float)a_p[j]*r[i-j];
		}
		
		k = -((float)r[i]+sum)/E[i-1];
		
		a_p[i] = k;
		
		for(uint8_t j=1; j<=i-1; j++)
		{
			a_n[j] = a_p[j] + k * a_p[i-j];
		}
		
		E[i] = (1.0-k*k) * E[i-1];
	}*/
	
	rv[0]=1.0;
	tmp[0]=1.0;
	rv[1]=-(float)r[1]/r[0];
	
	v=r[0]-(float)r[1]*r[1]/r[0];
	
	for(uint8_t j=1; j<=10; j++)
	{
		e=0.0;
		
		for(uint8_t i=0; i<=j-1; i++)
		{
			e+=r[j-i]*a[i-1];
			c=e/v;
			v=v-c*e;
			
		}
	}
	
	for(uint8_t i=0; i<11; i++)
		printf("%f\n", rv[i]);
	/*
	#ifdef DEBUG
	for(uint8_t i=0; i<10; i++)
	{
		for(uint8_t j=0; j<10; j++)
		{
			if(j<9)
				printf("%lld,", r[abs(i-j)]);
			else
				printf("%lld;\n", r[abs(i-j)]);
		}
	}
	
	for(uint8_t i=1; i<=10; i++)
			printf("%lld\n", r[i]);
	#endif
	*/
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
	
	LD_Solver(r, NULL);
	
	return 0;
}

