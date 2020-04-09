//Standard includes
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI			3.14159265358979323846
#endif

//#define DEBUG										//comment it out later
//define OVF_INFO
#define ERRORS

//-----------------------------ACELP includes & defines------------------------------
#include "LSP_codebooks.h"
//Global consts
#define FRAME_LEN		20							//voice frame length in ms
#define FRAME_SIZ		160							//voice frame samples number, 0.02s*8000Hz
#define	WINDOW_LEN		30							//length of window for autocorrelation computation
#define WINDOW_SIZ		240							//window samples number, 0.03s*8000Hz
#define ALPHA			(float)32735.0/32768.0		//alpha coeff for the pre-processing filter
#define GRID_SIZ		60							//grid granularity for LSP computation

//Global vars
float		w[WINDOW_SIZ];							//modified Hamming window w(n) coeffs for speech analysis
float		grid[GRID_SIZ];							//grid of values for LSP computation

uint64_t	frame = 0;								//frame number - for our info

//----------------------------------ACELP functions----------------------------------
//generate grid of values for LSP computation
void Grid_Generate(float *g)
{
	g[0] = 1.0;
	g[GRID_SIZ] = -1.0;
	
	for(uint8_t i=1; i<GRID_SIZ; i++)
		g[i] = cos((M_PI*i)/GRID_SIZ);
}

//LSP F(z) polynomial evaluation using Chebyshev polynomials
//arg1: input value, arg2: f() coeffs
//retval: evaluated value
float Chebyshev_Eval(float x, float *f)
{
	uint8_t n=5;	//coeffs num
	
	float b0, b1, b2;
	
	b2 = f[0];
	b1 = 2.0*x + f[1];
	
	for(uint8_t i=2; i<n; i++)
	{
		b0 = 2.0*x*b1 - b2 + f[i];
		b2 = b1;
		b1 = b0;
	}
	
	return x*b1 - b2 + 0.5*f[n];
}

//input speech signal pre-processing
//y[i] = x[i]/2 - x[i-1]/2 + alpha * y[i-1]
//arg1: present frame, arg2: output
void Speech_Pre_Process(int16_t *inp, int16_t *outp)
{
	#ifdef ERRORS
	if(inp==NULL || outp==NULL)
	{
		printf("\nNULL pointer at Speech_Pre_Process()\n");
		exit(0);
	}
	#endif
	
	outp[0]=inp[0]/2;
	
	for(uint8_t i=1; i<WINDOW_SIZ; i++)
	{
		outp[i] = inp[i]/2 - inp[i-1]/2 + ALPHA*outp[i-1];
	}
}

//compute the modified Hamming window w(n) coeffs for speech analysis
void Analysis_Window_Init(float *w)
{
	uint8_t L2 = 40;				//40 samples look ahead
	uint8_t L1 = WINDOW_SIZ-L2;		//200 samples
	
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
void Autocorr(int16_t *spch, int32_t *r)
{
	uint8_t ovf;
	int64_t tmp;				//for a[0] computation
	uint8_t norm_shift = 0;		//shifts left needed to normalize r[]

	//initially, set r[0]=1 to avoid r[] containing only zeros
	//and zero out the rest
	r[0]=1;
	memset(&r[1], 0, 10*sizeof(int32_t));
	
	//r[0] calculation
	//if r[0] overflows int32_t, divide the signal by 4
	do
	{
		ovf = 0;
		tmp = 0;
		
		for(uint8_t i=0; i<WINDOW_SIZ; i++)
		{
			tmp += (int64_t)spch[i] * (int64_t)spch[i];
		
			if(tmp > (int64_t)INT32_MAX)	//overflow occured?
			{
				//divide the signal by 4
				for(uint8_t j=0; j<WINDOW_SIZ; j++)
					spch[j] /= 4;
				
				ovf = 1;
				
				#ifdef OVF_INFO
				printf("Overflow occured in Autocorr() at i=%d\n", i);
				printf("val=%lld\n", tmp);
				#endif
				
				//break the "for" loop
				break;
			}
		}
	}
	while(ovf);
	
	r[0] = (int32_t)tmp;
	
	//r[0] normalization to the int32_t limit
	//multiply by 2 until we can't no more
	for(uint8_t i=0; i<32; i++)
	{
		while(r[0] < (INT32_MAX/2-1))
		{
			r[0]*=2;
			norm_shift++;
		}
	}
	
	//r[1]..r[10] calculation
	for(uint8_t i=1; i<=10; i++)
	{
		for(uint8_t j=0; j<WINDOW_SIZ; j++)
			r[i] += spch[j] * spch[j-i];
			
		//normalize
		for(uint8_t j=0; j<norm_shift; j++)
			r[i] *= 2;
	}
	
	#ifdef DEBUG
	printf("\nnorm_shift=%d\n", norm_shift);
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
void LD_Solver(int32_t *r, float *a)
{
	//previous frame coeffs - static variable!
	static float prev_a[11] = {1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	float k, alpha;
	float at[11], an[11];	//t: this iteration, n: next iteration
	float sum;
	
	k = -(float)r[1] / r[0];
	at[1] = k;
	alpha = (float)r[0] * (1.0-k*k);
	
	for(uint8_t i=2; i<=10; i++)
	{
		sum = 0.0;
		
		for(uint8_t j=1; j<=i-1; j++)
		{
			sum += (float)r[j]*at[i-j];
		}
		
		sum += (float)r[i];
		k = -sum / alpha;
		
		//test for filter stability
		//if case of instability, use previous coeffs
		//TODO: something might be wrong with this check
		if(fabs(k) > 32750.0/32767.0)	//close enough to 1.0
		{
			memcpy(a, prev_a, 11*sizeof(float));
			#ifdef ERRORS
			printf("Unstable filter, k=%1.2f at i=%d\n", k, i);
			#endif
			return;
		}
		
		//compute new coeffs
		for(uint8_t j=1; j<=i-1; j++)
		{
			an[j] = at[j] + k*at[i-j];
		}
		
		an[i] = k;
		alpha *= (1.0-k*k);
		
		memcpy(at, an, 11*sizeof(float));
	}
	
	prev_a[0]=1.0;
	memcpy(&prev_a[1], &at[1], 10*sizeof(float));
	
	//return solution
	if(a!=NULL)
	{
		a[0]=1.0;	//denominator of the H(z) is A(z)=1+sum(a*z)
		memcpy(&a[1], &at[1], 10*sizeof(float));
	}
	
	#ifdef DEBUG
	for(uint8_t i=0; i<11; i++)
		printf("a[%d]=%f\n", i, a[i]);
	#endif
}

//convert LP to LSP
//arg1: previous frame LSP array, arg2: present frame LP array, arg3: output LSP array
void LP_LSP(float *prev_LSP, float *a, float *LSP)
{
	float f1[6] = {1.0, 0, 0, 0, 0, 0};
	float f2[6] = {1.0, 0, 0, 0, 0, 0};
	float *coefs;							//coeff set that we are using
	uint8_t found=0;						//found roots
	uint8_t loc=0;							//location in the grid
	float x1, x2, y1, y2, xm, ym, x, y;		//vals for root search

	//5 polynomial coeffs
	for(uint8_t i=0; i<5; i++)
 	{
		f1[i+1] = a[i+1] + a[10-i] - f1[i];
		f2[i+1] = a[i+1] - a[10-i] + f2[i];
	}
	
	#ifdef DEBUG
	//evaluation of the polynomial - root search
	printf("\n");
	for(uint8_t i=0; i<GRID_SIZ; i++)
		;//printf("%f\n", Chebyshev_Eval(grid[i], f2));
	#endif
	
	//look for roots in f1 first
	coefs = f1;
	
	//search init
	x1 = grid[0];
	y1 = Chebyshev_Eval(x1, coefs);
	
	//search for the roots
	//until we have 10 or we have searched thru all the grid values (0..pi)
	while(found<10 && loc<GRID_SIZ)
	{
		loc++;	//move thru the grid
		
		x2 = x1;
		y2 = y1;
		x1 = grid[loc];
		y1 = Chebyshev_Eval(x1, coefs);
		
		//check for a sign change
		if(y1*y2 <= 0)
		{
			//divide the range 4 times
			for(uint8_t i=0; i<4; i++)
      		{
        		xm = 0.5 * (x1+x2);
        		ym = Chebyshev_Eval(xm, coefs);
  				
  				//sign change in the lower half?
				if(y1*ym <= 0)
        		{
        			//move there
          			y2 = ym;
          			x2 = xm;
        		}
        		//same thing here - zero crossing in the second half?
        		else
        		{
        			//move there
          			y1 = ym;
          			x1 = xm;
        		}
        	}
        	
        	//linear interpolation for the fine root value
        	x = x2-x1;
        	y = y2-y1;
        	
			if(fabs(y)<0.0001)	//unsafe to compare floats to 0.0
				x = x1;
			else
			{
				y = (x2-x1)/(y2-y1);
				y=fabs(y);
				x1 = x1 - y1*y;
			}
			
			LSP[found]=x1;
			found++;
			
			#ifdef DEBUG
			printf("%f|%f|%d|%f\n", x1, y2, loc, 8000.0/(2*M_PI)*acos(x1));
			#endif
			
			//swap f1 with f2 and vice-versa, for next search
			if(coefs == f1)
			{
				coefs = f2;
			}
			else
			{
				coefs = f1;
  	    	}
		}
        	
		//apply new value of y1
		y1 = Chebyshev_Eval(x1, coefs);
	}
	
	//check if we have found all 10 roots
	//if not - copy old roots and use them
	if(found<10 && prev_LSP!=NULL)
	{
		memcpy(LSP, prev_LSP, 10*sizeof(float));
		#ifdef DEBUG
		printf("\nLess than 10 roots found in LP_LSP()\n");
		#endif
	}
}

//split vector quantization of LSP parameters
//full codebook search with squared error metric (saving one division)
//arg1: LSPs in cosine domain (10), arg2: quantized LSPs output (10), arg3: codebook indices output (3)
void LSP_SVQ(float *lsp, float *q_lsp, uint16_t *ind)
{
	uint16_t ind_rv[3];
	
	float se;
	float delta;
	float min=10000.0;
	
	//codebook 1 search
	for(uint16_t i=0; i<size_cb1; i++)
	{
		se=0.0;
		
		for(uint8_t j=0; j<3; j++)
		{
			delta = lsp[j]-cb1[i*3+j];
			se += delta*delta;
		}
		
		if(se < min)
		{
			min = se;
			ind_rv[0]=i;
		}
	}
	    
	min   =10000.0;
	
	//codebook 2 search
	for(uint16_t i=0; i<size_cb2; i++)
	{		
		se=0.0;
		
		for(uint8_t j=0; j<3; j++)
		{
			delta = lsp[3+j]-cb2[i*3+j];
			se += delta*delta;
		}
		
		if(se < min)
		{
			min = se;
			ind_rv[1]=i;
		}
	}
	  
	min   =10000.0;
	
	//codebook 3 search
	for(uint16_t i=0; i<size_cb3; i++)
	{
		se=0.0;
		
		for(uint8_t j=0; j<4; j++)
		{
			delta = lsp[6+j]-cb3[i*4+j];
			se += delta*delta;
		}
		
		if(se < min)
		{
			min = se;
			ind_rv[2]=i;
		}
	}
	
	//return quantized vector...
	memcpy(&q_lsp[0], &cb1[ind_rv[0]], 3*sizeof(float));
	memcpy(&q_lsp[3], &cb2[ind_rv[1]], 3*sizeof(float));
	memcpy(&q_lsp[6], &cb3[ind_rv[2]], 4*sizeof(float));
	
	//...and codebook indices
	memcpy(ind, ind_rv, 3*sizeof(uint16_t));
}

//convert LSP coeffs to F1(z) or F2(z)
//arg1: LSP array of length 10, arg2: F_1(z) or F_2(z) coefficients output
void LSP_Poly(float *lsp, float *f)
{
	uint8_t k=0;
	
	f[0] = 1.0;
	f[1] = -2.0 * lsp[k];
	k+=2;
	
	for(uint8_t i=2; i<=5; i++)
	{
		f[i] = -2.0*lsp[k]*f[i-1]+2.0*f[i-2];
		
		for(int8_t j=i-1; j>=1; j--)
		{
			if(j>1)
				f[j] = f[j] - 2.0*lsp[k]*f[j-1] + f[j-2];
			else
				f[j] = f[j] - 2.0*lsp[k]*f[j-1];	//f(-1)=0
		}
		
		k+=2;
	}
}

//convert LSP to LP (A(z))
//arg1: input LSP array, arg2: computed LP array
void LSP_LP(float *lsp, float *a)
{
	float f1[6], f2[6];
	
	//get F1(z) and F2(z) coeffs
	LSP_Poly(&lsp[0], f1);
	LSP_Poly(&lsp[1], f2);
	
	for(int8_t i=5; i>0; i--)
	{
		f1[i] += f1[i-1];
		f2[i] -= f2[i-1];
	}
	
	a[0]=1.0;
	
	int8_t i, j;
	for(i=1, j=10; i<=5; i++, j--)
	{
		a[i] = 0.5 * (f1[i] + f2[i]);
		a[j] = 0.5 * (f1[i] - f2[i]);
	}
}

//encode voice frame
//arg1: input speech samples, 16-bit signed integer, arg2: 137 unpacked output bits
void ACELP_EncodeFrame(int16_t *speech, uint8_t *out)
{
	//first call of this function?
	//make it global later,
	//so it can be accessed outside of this function
	static uint8_t first=1;
	
	//local buffers for speech frame manipulation
	int16_t spch_in[WINDOW_SIZ];
	int16_t spch_out[WINDOW_SIZ];
	
	int32_t		r[11];									//autocorrelation values
	float		lp[11];									//LP coeffs (10, but starting from lp[1], lp[0]=1.0)
	float		lsp[10];								//LSP coeffs in cosine domain
	
	//quantized and unquantized LSP vectors from the previous frame and this one
	static float q_lsp_prev[10];
	float q_lsp_this[10];
	static float lsp_prev[10];
	float lsp_this[10];
	
	//LSP codebook indices for this frame
	uint16_t lsp_cb_indices[3];
	
	//quantized LSP vector interpolation for 4 subframes
	float q_lsp_1[10];
	float q_lsp_2[10];
	float q_lsp_3[10];
	float q_lsp_4[10];
	
	//"real", unquantized LSP vector interpolation for 4 subframes
	float lsp_1[10];
	float lsp_2[10];
	float lsp_3[10];
	float lsp_4[10];
	
	//first frame? initialize "previous" frame LSP vectors
	if(first)
	{
		//previous quantized LSP vector
		lsp_prev[0] = q_lsp_prev[0] = 30000.0/32768.0;
		lsp_prev[1] = q_lsp_prev[1] = 26000.0/32768.0;
		lsp_prev[2] = q_lsp_prev[2] = 21000.0/32768.0;
		lsp_prev[3] = q_lsp_prev[3] = 15000.0/32768.0;
		lsp_prev[4] = q_lsp_prev[4] = 8000.0/32768.0;
		lsp_prev[5] = q_lsp_prev[5] = 0.0;
		lsp_prev[6] = q_lsp_prev[6] = -8000.0/32768.0;
		lsp_prev[7] = q_lsp_prev[7] = -15000.0/32768.0;
		lsp_prev[8] = q_lsp_prev[8] = -21000.0/32768.0;
		lsp_prev[9] = q_lsp_prev[9] = -26000.0/32768.0;
		
		first = 0;
	}
	
	//pre processing and windowing
	Speech_Pre_Process(speech, spch_out);
	memcpy(spch_in, spch_out, WINDOW_SIZ*sizeof(int16_t));
	Window_Speech(spch_in, spch_out);
	
	//compute LSPs
	Autocorr(spch_out, r);
	LD_Solver(r, lp);
	LP_LSP(lsp_prev, lp, lsp);
	
	//LSPs now
	memcpy(lsp_this, lsp, 10*sizeof(float));
	
	//quantize LSPs
	LSP_SVQ(lsp, q_lsp_this, lsp_cb_indices);
	
	//interpolate quantized LSP vector for subframes 4,3,2,1
	memcpy(q_lsp_4, q_lsp_this, 10*sizeof(float));
	for(uint8_t i=0; i<10; i++)
	{
		q_lsp_3[i] = 0.75*q_lsp_4[i] + 0.25*q_lsp_prev[i];
		q_lsp_2[i] = 0.50*q_lsp_4[i] + 0.50*q_lsp_prev[i];
		q_lsp_1[i] = 0.25*q_lsp_4[i] + 0.75*q_lsp_prev[i];
	}
	
	//interpolate unquantized LSP vector for subframes 4,3,2,1
	memcpy(lsp_4, lsp_this, 10*sizeof(float));
	for(uint8_t i=0; i<10; i++)
	{
		lsp_3[i] = 0.75*lsp_4[i] + 0.25*lsp_prev[i];
		lsp_2[i] = 0.50*lsp_4[i] + 0.50*lsp_prev[i];
		lsp_1[i] = 0.25*lsp_4[i] + 0.75*lsp_prev[i];
	}
	
	//now we have both quantized and unquantized LSP vectors for further computations
	//we can change them back to {a_i} for the A(z)
	
	//LSP to A(z) conversion
	float a[11];
	LSP_LP(lsp_1, a);
	
	
	
	//update LSPs
	memcpy(q_lsp_prev, q_lsp_this, 10*sizeof(float));
	memcpy(  lsp_prev,   lsp_this, 10*sizeof(float));
}

//initialize consts
void ACELP_Init(float *search_grid, float *analysis_window)
{
	Grid_Generate(search_grid);
	Analysis_Window_Init(analysis_window);
}

//main routine
//argv[1]: file name (RAW, signed 16-bit, Little-Endian, 8000Hz)
int main(uint8_t argc, uint8_t *argv[])
{
	FILE *aud;
	int16_t spch[WINDOW_SIZ];
	
	if(argc==2)
	{
		printf("Loading \"%s\"...\n\n", argv[1]);
		
		aud = fopen(argv[1], "rb");
		
		if(aud==NULL)
		{
			printf("No file named \"%s\"\n", argv[1]);
			printf("Exiting.\n");
			return 1;
		}
		else
		{
			//initialize consts etc.
			ACELP_Init(grid, w);
			
			//load 30ms frames, overlapping
			while(fread(spch, 2, WINDOW_SIZ, aud)==WINDOW_SIZ)
			{
				frame++;
				
				//take us 10ms back (80 samples * sizeof(int16_t))
				fseek(aud, -160, 1);
				
				printf("Frame %d\n", frame);
				
				ACELP_EncodeFrame(spch, NULL);
			}
		}
	}
	else
	{
		printf("Invalid params\nExiting.\n");
		return 1;
	}
	
	return 0;
}

