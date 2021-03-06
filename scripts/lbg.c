#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//enough to hold a 4-D vector
struct vec
{
	float x;
	float y;
	float z;
	float a;
};

//input data set
struct vec data_in[10000000];

//codebooks (arrays of vectors)
struct vec q1[256];
struct vec q2[512];
struct vec q3[512];

//512 dynamic length arrays for Voronoi cell vector binning ("partitions")
typedef struct
{
	struct vec *array;
	size_t used;
	size_t size;
} Array;

Array partition[512];
uint32_t part_count[512];		//their element counters

void initArray(Array *a, size_t initialSize)
{
	a->array = (struct vec*)malloc(initialSize*sizeof(struct vec));
	a->used = 0;
	a->size = initialSize;
}

void insertArray(Array *a, struct vec *element)
{
	if(a->used == a->size)
	{
		a->size++;
		a->array = (struct vec*)realloc(a->array, a->size*sizeof(struct vec));
	}
	
	a->array[a->used++] = *element;
}

void freeArray(Array *a)
{
	free(a->array);
	a->array = NULL;
	a->used = a->size = 0;
}
//------------------------------------------------------------------------

//distance (squared) between vectors
//arg1: input vector, arg2: input vector, arg3: dimension = {3, 4}
float distance(struct vec *inp_vector1, struct vec *inp_vector2, uint8_t dimension)
{
	if(dimension==3)
	{
		float x2=inp_vector2->x;
		float y2=inp_vector2->y;
		float z2=inp_vector2->z;
		
		float x1=inp_vector1->x;
		float y1=inp_vector1->y;
		float z1=inp_vector1->z;
		
		return ((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	}
	else if(dimension==4)
	{
		float x2=inp_vector2->x;
		float y2=inp_vector2->y;
		float z2=inp_vector2->z;
		float a2=inp_vector2->a;
		
		float x1=inp_vector1->x;
		float y1=inp_vector1->y;
		float z1=inp_vector1->z;
		float a1=inp_vector1->a;
		
		return ((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)+(a2-a1)*(a2-a1));
	}
}

//average distance (squared) between a set of vectors and a vector
//arg1: input vector, arg2: input vector set, arg3: input vector set count, arg4: dimension = {3, 4}
float distance_set(struct vec *inp_vector, struct vec *inp_data, uint32_t count, uint8_t dimension)
{
	float tmp=0.0;
	
	for(uint32_t i=0; i<count; i++)
	{
		tmp += distance(inp_vector, &inp_data[i], dimension);
		//printf("tmp=%f\n", tmp);
	}
		
	return tmp/count;
}

//compute centroid of a set of vectors
//arg1: output vector, arg2: vector count, arg3: set of vectors, arg4: dimension
void centroid(struct vec *out, struct vec *inp_data, uint32_t count, uint8_t dimension)
{
	float sum_x=0.0;
	float sum_y=0.0;
	float sum_z=0.0;
	float sum_a=0.0;
	
	if(dimension==3)
	{
		for(uint32_t i=0; i<count; i++)
		{
			sum_x += inp_data[i].x;
			sum_y += inp_data[i].y;
			sum_z += inp_data[i].z;
		}
		
		out->x = sum_x/count;
		out->y = sum_y/count;
		out->z = sum_z/count;
	}
	else
	{
		for(uint32_t i=0; i<count; i++)
		{
			sum_x += inp_data[i].x;
			sum_y += inp_data[i].y;
			sum_z += inp_data[i].z;
			sum_a += inp_data[i].a;
		}
		
		out->x = sum_x/count;
		out->y = sum_y/count;
		out->z = sum_z/count;
		out->a = sum_a/count;
	}
}

//multiply vector by (1.0+e)
void vec_mult(struct vec *out, float e, struct vec *inp)
{	
	out->x = inp->x * (1.0 + e);
	out->y = inp->y * (1.0 + e);
	out->z = inp->z * (1.0 + e);
	out->a = inp->a * (1.0 + e);
}

int main(void)
{
	uint32_t data_count;
	uint16_t cell_count=1;
	float epsilon = 0.0001;
	
	memset(part_count, 0, 512*sizeof(uint32_t));

	//test data set
	data_in[0].x=0.99;
	data_in[0].y=0.0;
	data_in[0].z=0.0;
	
	data_in[1].x=0.0;
	data_in[1].y=1.01;
	data_in[1].z=0.0;
	
	data_count=2;
	
	//codebook ^q1 = {q1, q2, q3}
	//initial iteration
	centroid(&q1[0], data_in, data_count, 3);
	
	printf("q[0] = ( %f, %f, %f )\n\n", q1[0].x, q1[0].y, q1[0].z);
	
	//first split
	vec_mult(&q1[1], epsilon, &q1[0]);
	vec_mult(&q1[0], -epsilon, &q1[0]);
	cell_count++;
	
	printf("q1[0] = ( %f, %f, %f )\n", q1[0].x, q1[0].y, q1[0].z);
	printf("q1[1] = ( %f, %f, %f )\n\n", q1[1].x, q1[1].y, q1[1].z);
	
	//iteration 2
	for(uint32_t i=0; i<data_count; i++)
	{
		float min=20e6;
		uint32_t c=0;	//which codeword is the closest
		
		//for every codeword
		for(uint16_t j=0; j<cell_count; j++)
		{		
			float tmp=distance(&q1[j], &data_in[i], 3);
				
			if(tmp < min)
			{
				min=tmp;
				c=j;
			}
		}
		
		//move input data to a "bin"
		insertArray(&partition[c], &data_in[i]);
	}
	
	//recalculate codewords
	centroid(&q1[2], &partition[0].array[0], partition[0].used, 3);
	centroid(&q1[3], &partition[1].array[0], partition[1].used, 3);
	
	vec_mult(&q1[1], epsilon, &q1[3]);
	vec_mult(&q1[0], -epsilon, &q1[2]);
	
	cell_count++;
	
	printf("q1[0] = ( %f, %f, %f )\n", q1[0].x, q1[0].y, q1[0].z);
	printf("q1[1] = ( %f, %f, %f )\n", q1[1].x, q1[1].y, q1[1].z);
	printf("q1[2] = ( %f, %f, %f )\n", q1[2].x, q1[2].y, q1[2].z);
	printf("q1[3] = ( %f, %f, %f )\n", q1[3].x, q1[3].y, q1[3].z);
	
	return 0;
}

