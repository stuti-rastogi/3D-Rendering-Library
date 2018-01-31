/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image = NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (reset) {          /* open and load texture file */
		fd = fopen("texture", "rb");
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (image == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		reset = 0;          /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */
	if (u < 0)
		u = 0;
	else if (u > 1)
		u = 1;

	if (v < 0)
		v = 0;
	else if (v > 1)
		v = 1;

	float u_scale = u * (xs - 1);
	float v_scale = v * (ys - 1);

	int u_low = floor(u_scale);
	int v_low = floor(v_scale);
	int u_high = ceil(u_scale);
	int v_high = ceil(v_scale);

	float s = u_scale - u_low;
	float t = v_scale - v_low;

	// A: (u_low, v_low)
	// B: (u_high, v_low)
	// C: (u_high, v_high)
	// D: (u_low, v_high)
	//Color(p) = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A 
	for (int i = 0; i < 3; i++)
	{
		color[i] = (s * t * image[v_high*xs + u_high][i]) +
			((1 - s) * t * image[v_high*xs + u_low][i]) +
			((1 - t) * s * image[v_low*xs + u_high][i]) +
			((1 - s) * (1 - t) * image[v_low*xs + u_low][i]);
		//color[i] = image[v_low*(xs) + u_low][i];
	}

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	/*Set X = [(u - W / 2) / (W / 2);
	(v - H / 2) / (H / 2)]*/

	float w = 512.0;
	float h = 512.0;

	float x_real = u;//(u - (w / 2)) / (w / 2);
	float x_im = v;// (v - (h / 2)) / (h / 2);

	float c_real = 0.45;
	float c_im = -0.527;
	float N = 200;
	float length = 0;
	int i;

	for (i = 0; i < N; i++)
	{
		x_real = pow(x_real, 2) - pow(x_im, 2) + c_real;
		x_im = (2 * x_real * x_im) + c_im;;
		length = sqrt(pow(x_real, 2) + pow(x_im, 2));
		if (length >= 2)
			break;
	}

	float z = i / N;

	//int N = 300;
	//float c_real = 0.123;//-0.12375;
	//float c_im = -0.015;//0.56805;
	//float x_real = u;
	//float x_im = v;
	//float length = 0;

	//for (int i = 0; i < N; i++)
	//{
	//	// x = u+vi
	//	// compute x^2
	//	x_real = pow(x_real,3) - (3*x_real*pow(x_im, 2)) + c_real;
	//	x_im = (3 * pow(x_real,2) * x_im) - pow(x_im,3) + c_im;

	//	/*x_real = pow(x_real,2) - pow(x_im,2) + c_real;
	//	x_im = (2*x_real*x_im) + c_im;*/
	//	length = sqrt(pow(x_real, 2) + pow(x_im, 2));
	//	if (length >= 2)
	//		break;
	//}

	// LUT
	float divisions[11] = { 0.0, 0.1, 0.2,
		0.3, 0.4, 0.5,
		0.6, 0.7, 0.8,
		0.9, 1.0 };
	GzColor colors[11] = { { 1.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0 },{ 0.980, 0.972, 0.239 },
	{ 0.098, 0.090, 0.529 },{ 0.576, 0.047, 0.654 },{ 0.756, 0.870, 0.129 },
	{ 1, 0.780, 0.019 },{ 0.019, 0.98, 0.968 },{ 0.074, 0.168, 0.486 },
	{ 0.541, 0.545, 0.517 },{ 0.07, 0.04, 0.46 } };

	int low = 0;
	int high = 1;
	if (z >= 0.0 && z < 0.1)
	{
		low = 0;
		high = 1;
	}
	else if (z >= 0.1 && z < 0.2)
	{
		low = 1;
		high = 2;
	}
	else if (z >= 0.2 && z < 0.3)
	{
		low = 2;
		high = 3;
	}
	else if (z >= 0.3 && z < 0.4)
	{
		low = 3;
		high = 4;
	}
	else if (z >= 0.4 && z < 0.5)
	{
		low = 4;
		high = 5;
	}
	else if (z >= 0.5 && z < 0.6)
	{
		low = 5;
		high = 6;
	}
	else if (z >= 0.6 && z < 0.7)
	{
		low = 6;
		high = 7;
	}
	else if (z >= 0.7 && z < 0.8)
	{
		low = 7;
		high = 8;
	}
	else if (z >= 0.8 && z < 0.9)
	{
		low = 8;
		high = 9;
	}
	else if (z >= 0.9 && z <= 1.0)
	{
		low = 9;
		high = 10;
	}

	float s_i = divisions[low];
	float s_i1 = divisions[high];

	float a = (s_i1 - z) / (s_i1 - s_i);
	float b = (z - s_i) / (s_i1 - s_i);

	for (int i = 0; i < 3; i++)
	{
		color[i] = (a * colors[low][i]) + (b * colors[high][i]);
	}

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if (image != NULL)
		free(image);
	return GZ_SUCCESS;
}

