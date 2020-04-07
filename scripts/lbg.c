#include <stdio.h>
#include <inttypes.h>
#include <math.h>

//enough to hold a 4-D vector
struct vec
{
	float x;
	float y;
	float z;
	float a;
};

struct vec data_in[10000000];						//input data set
struct vec *partition[512];							//pointers to binned vectors (for Voronoi cell vector binning)

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
		tmp += distance(inp_vector, &inp_data[i], dimension)/count;
	}
		
	return tmp;
}

//compute centroid of a set of vectors
//arg1: set of vectors, arg2: vector count, arg3: output vector, arg4: dimension
void centroid(struct vec *inp_data, uint32_t count, struct vec *out, uint8_t dimension)
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
void vec_mult(struct vec *inp, float e)
{
	inp->x *= (1.0 + e);
	inp->y *= (1.0 + e);
	inp->z *= (1.0 + e);
	inp->a *= (1.0 + e);
}

int main(void)
{
	data_in[0].x=1.0;
	data_in[0].y=0.0;
	data_in[0].z=0.0;
	
	data_in[1].x=2.0;
	data_in[1].y=2.0;
	data_in[1].z=2.0;
	
	struct vec test_vec={1.0, 1.0, 1.0, 0.0};

	printf("%f\n", distance_set(&test_vec, data_in, 2, 3));
	
	return 0;
}

