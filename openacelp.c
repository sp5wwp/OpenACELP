//Standard includes
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

//-----------------------------ACELP includes & defines------------------------------
#define MAX_32	 2147483647
#define MIN_32	-2147483648
#define HALF_32	 1073741824

#define MAX_16	 32767
#define MIN_16	-32768
#define HALF_16	 16384


//Global consts
#include "N_POW2.h"	//powers of 2 array
#include "LOG2.h"	//array for Log2 computation	

//Global flags
uint8_t ovf=0;		//global overflow flag
uint8_t carry=0;	//global carry flag

//----------------------------------Basic functions----------------------------------
//---------------------------------16-bit arithmetic---------------------------------
//cast an int32_t into int16_t with saturation
int16_t saturate(int32_t val)
{
	int16_t rv;

	if(val > MAX_16)
	{
		ovf = 1;
		return MAX_16;
	}
	else if(val < MIN_16)
	{
		ovf = 1;
		return MIN_16;
	}
	else
	{
		ovf = 0;
		return (int16_t)val;
	}
}

//absolute value with saturation
int16_t abs_s(int16_t val)
{
	if(val==MIN_16)
		return MAX_16;
	else
		return (int16_t)abs(val);
}

//sum with saturation
int16_t add_s(int16_t a, int16_t b)
{
	return saturate((int32_t)a + (int32_t)b);
}

//subtraction with saturation
int16_t sub(int16_t a, int16_t b)
{
	return saturate((int32_t)a - (int32_t)b);
}

//calculate val1*val2, scaled
int16_t mult(int16_t val1, int16_t val2)
{
	return saturate(((int32_t)val1 * (int32_t)val2) / (1<<15));
}

//---------------------------------32-bit arithmetic---------------------------------
int32_t L_shl(int32_t val1, int16_t val2);
int32_t L_shr(int32_t val1, int16_t val2);

//sum with saturation
int32_t L_add_s(int32_t val1, int32_t val2)
{
	int64_t rv = val1 + val2;

	if(rv > MAX_32)
	{
		ovf = 1;
		return MAX_32;
	}
	else if(rv < MIN_32)
	{
		ovf = 1;
		return MIN_32;
	}
	else
	{
		return (int32_t)rv;
	}
}

//subtraction with saturation (val1-val2)
int32_t L_sub_s(int32_t val1, int32_t val2)
{
	int64_t rv = (int64_t)val1 - (int64_t)val2;

	if(rv > MAX_32)
	{
		ovf = 1;
		return MAX_32;
	}
	else if(rv < MIN_32)
	{
		ovf = 1;
		return MIN_32;
	}
	else
	{
		return (int32_t)rv;
	}
}

//L_mult is the 32 bit result of the multiplication of var1 times var2
int32_t L_mult(int16_t val1, int16_t val2)
{
	int64_t rv = ((int32_t)val1 * (int32_t)val2);

	if(rv > MAX_32)
	{
		ovf = 1;
		return MAX_32;
	}
	else if(rv < MIN_32)
	{
		ovf = 1;
		return MIN_32;
	}
	else
		return (int32_t)rv;
}

//L_mult_shl is the 32 bit result of the multiplication of var1 times var2 with one shift left
int32_t L_mult_shl(int16_t val1, int16_t val2)
{
	int64_t rv = ((int64_t)val1 * (int64_t)val2)*2;

	if(rv > MAX_32)
	{
		ovf = 1;
		return MAX_32;
	}
	else if(rv < MIN_32)
	{
		ovf = 1;
		return MIN_32;
	}
	else
		return (int32_t)rv;
}

//calculate val1*val2+val3 (multiply-accumulate)
int32_t L_mac(int16_t val1, int16_t val2, int32_t val3)
{
	return L_add_s(L_mult(val1, val2), val3);
}

//calculate shl(val1*val2, 1)+val3 (multiply-accumulate with shift left)
int32_t L_mac_shl(int16_t val1, int16_t val2, int32_t val3)
{
	return L_add_s(L_mult_shl(val1, val2), val3);
}

//arithmetically shift the 32 bit input val1 left by val2 positions, with saturation
int32_t L_shl(int32_t val1, int16_t val2)
{
	if(val2 <= 0)
		return L_shr(val1, -val2);
	else
	{
		int64_t rv = (int64_t)val1 * (1<<(uint16_t)val2);
		
		if(rv > MAX_32)
		{
			ovf = 1;
			return MAX_32;
		}
		else if(rv < MIN_32)
		{
			ovf = 1;
			return MIN_32;
		}
		else
			return (int32_t)rv;
	}
}

//arithmetically shift the 32 bit input val1 right by val2 positions, with saturation
int32_t L_shr(int32_t val1, int16_t val2)
{
	if(val2 < 0)
		return L_shl(val1, -val2);
	else
	{
		int64_t rv = (int64_t)val1 / (1<<(uint16_t)val2);
		
		if(rv > MAX_32)
		{
			ovf = 1;
			return MAX_32;
		}
		else if(rv < MIN_32)
		{
			ovf = 1;
			return MIN_32;
		}
		else
			return (int32_t)rv;
	}
}

//calculate shl(val1*val2, 1)-val3 with saturation
int32_t L_msu_shl(int16_t val1, int16_t val2, int32_t val3)
{
	return L_sub_s(val3, L_mult_shl(val1, val2));
}

//calculate val1*val2-val3 with saturation
int32_t L_msu(int16_t val1, int16_t val2, int32_t val3)
{
	return L_sub_s(val3, L_mult(val1, val2));
}

//"extract the lower portion of int32_t"
//TODO: this function might be magic, investigate its behaviour. it just casts int32_t into int16_t
int16_t extract_l(int32_t val)
{
	return (int16_t)val;
}

//calculate shl(val1, shift)+val2
int32_t add_shl(int16_t val1, int16_t shift, int32_t val2)
{
	return L_msu(val1, N_POW2[shift], val2);
}
 
//calculate shl(val1, 16)+val2
int32_t add_shl16(int16_t val1, int32_t val2)
{
	return L_msu_shl(val1, -32768, val2);
}

//calculate shl(val1, shift)-val2
int32_t sub_shl(int16_t val1, int16_t shift, int32_t val2)
{
	return L_mac(val1, N_POW2[shift], val2);
}

//calculate shl(val1, 16)-val2
int32_t sub_shl16(int16_t val1, int16_t shift, int32_t val2)
{
	return L_mac_shl(val1, -32768, val2);
}

//calculate -val with saturation
int16_t negate(int16_t val)
{
	if(val==MIN_16)
		return MAX_16;
	else
		return -val;
}

//bit shifts needed to normalize int32_t to range MIN_32..-HALF_32 and HALF_32..MAX_32
int16_t norm_l(int32_t val)
{
	if(val == 0)
	{
		return 0;
	}
	else
	{
		if(val == -1)
		{
			return 31;
		}
		else
		{
			if(val < 0)
				val = ~val;

			int16_t rv=0;

			for(; val < HALF_32; rv++)
			{
				val <<= 1;
			}
			
			return rv;
		}
	}
}

//bit shifts needed to normalize int16_t to range MIN_16..-HALF_16 and HALF_16..MAX_16
int16_t norm_s(int16_t val)
{
	if(val == 0)
	{
		return 0;
	}
	else
	{
		if(val == -1)
		{
			return 15;
		}
		else
		{
			if(val < 0)
				val = ~val;

			int16_t rv=0;

			for(; val < HALF_16; rv++)
			{
				val <<= 1;
			}
			
			return rv;
		}
	}
}

//variable normalisation of a 32 bit integer (val1)
//val2 gives the maximum number of left shift to be done
//val3 returns the actual number of left shift
int32_t norm_v(int32_t val1, int16_t val2, int16_t *val3)
{
	int16_t shift = norm_l(val1);

	if(sub(shift, val2) > 0)
		shift = val2;
	*val3 = shift;

	return L_shl(val1, shift);
}

//store high part of val1 with a left shift of val2
int16_t store_h(int32_t val1, int16_t val2)
{
	static const int16_t SHR[8] = {16, 15, 14, 13, 12, 11, 10, 9};

	return extract_l(L_shr(val1, SHR[val2]));
}

//load the 16 bit val with a left shift of 16 into 32 bit output
int32_t Load_shl16(int16_t val)
{
	return L_msu_shl(val, -32768, 0);
}

//----------------------------------Mid-level functions----------------------------------
//innovative codebook search
int16_t Innov_Codebook_Search(int16_t dn[], int16_t f[], int16_t h[], int16_t rr[][32],
             int16_t cod[], int16_t y[], int16_t *sign, int16_t *shift_code)
{
	;
}

//LPC coefficients calculation
//based on the Levison-Durbin algorithm
//in double precision
void LD_solver_32(int16_t Rh[], int16_t Rl[], int16_t A[])
{
	;
}

int main(void)
{
	int16_t val=norm_l(-2);
	
	printf("%ld, ovf=%d\n", val, ovf);
	
	return 0;
}

