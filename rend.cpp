/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180.0;
	GzMatrix result =
	{
		1, 0, 0, 0,
		0, cos(radians), -sin(radians), 0,
		0, sin(radians), cos(radians), 0,
		0, 0, 0, 1
	};

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = result[i][j];
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180.0;
	GzMatrix result =
	{
		cos(radians), 0, sin(radians), 0,
		0, 1, 0, 0,
		-sin(radians), 0, cos(radians), 0,
		0, 0, 0, 1
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = result[i][j];
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180.0;
	GzMatrix result =
	{
		cos(radians), -sin(radians), 0, 0,
		sin(radians), cos(radians), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = result[i][j];
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	GzMatrix result =
	{
		1, 0, 0, translate[0],
		0, 1, 0, translate[1],
		0, 0, 1, translate[2],
		0, 0, 0, 1
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = result[i][j];
		}
	}
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	GzMatrix result =
	{
		scale[0], 0, 0, 0,
		0, scale[1], 0, 0,
		0, 0, scale[2], 0,
		0, 0, 0, 1
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = result[i][j];
		}
	}
	return GZ_SUCCESS;
}

GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	-- set display resolution
	-- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	-- allocate memory for pixel buffer
	*/
	xres = xRes;
	yres = yRes;
	int resolution = xres * yres;
	int framebufferSize = 3 * resolution;
	framebuffer = (char*)malloc(3 * sizeof(char) * xRes * yRes);
	pixelbuffer = new GzPixel[resolution];

	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/
	matlevel = 0;
	//default camera
	GzCamera defaultCamera;
	defaultCamera.position[X] = DEFAULT_IM_X;
	defaultCamera.position[Y] = DEFAULT_IM_Y;
	defaultCamera.position[Z] = DEFAULT_IM_Z;

	defaultCamera.lookat[X] = 0.0;
	defaultCamera.lookat[Y] = 0.0;
	defaultCamera.lookat[Z] = 0.0;

	defaultCamera.worldup[X] = 0.0;
	defaultCamera.worldup[Y] = 1.0;
	defaultCamera.worldup[Z] = 0.0;

	defaultCamera.FOV = DEFAULT_FOV;
	GzPutCamera(defaultCamera);

	//Xsp
	GzMatrix xsp_temp =
	{
		(xres / 2), 0, 0, (xres / 2),
		0, (-yres / 2), 0, (yres / 2),
		0, 0, MAXINT, 0,
		0, 0, 0, 1
	};

	memmove(Xsp, xsp_temp, sizeof(Xsp));
}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	delete[] framebuffer;
	delete[] pixelbuffer;
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	for (int i = 0; i < xres; i++)
	{
		for (int j = 0; j < yres; j++)
		{
			GzPut(i, j, 1991, 1799, 1574, 1, MAXINT);		// some default value
			//GzPut(i, j, 4095, 4095, 4095, 1, MAXINT);
		}
	}
	return GZ_SUCCESS;
}

float GzRender::calculateNorm(GzCoord vector)
{
	float norm, x, y, z;
	x = vector[0];
	y = vector[1];
	z = vector[2];
	norm = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	return norm;
}

void GzRender::calculateCrossProduct(GzCoord a, GzCoord b, float *resultX, float *resultY, float *resultZ)
{
	*resultX = (a[1] * b[2]) - (a[2] * b[1]);
	*resultY = (a[2] * b[0]) - (a[0] * b[2]);
	*resultZ = (a[0] * b[1]) - (a[1] * b[0]);
}

float GzRender::calculateDotProduct(GzCoord a, GzCoord b)
{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void GzRender::scalarVectorMultiply(float a, GzCoord v, float *resultX, float *resultY, float *resultZ)
{
	*resultX = a * v[0];
	*resultY = a * v[1];
	*resultZ = a * v[2];
}

void GzRender::calculateDifferenceVectors(GzCoord a, GzCoord b, float *resultX, float *resultY, float *resultZ)
{
	*resultX = a[0] - b[0];
	*resultY = a[1] - b[1];
	*resultZ = a[2] - b[2];
}

void GzRender::multiplyMatrices(GzMatrix mat1, GzMatrix mat2, GzMatrix *product)
{
	GzMatrix result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			float sum = 0;
			for (int k = 0; k < 4; k++)
			{
				sum = sum + (mat1[i][k] * mat2[k][j]);
			}
			result[i][j] = sum;
		}
	}
	memmove(*product, result, sizeof(result));
}

int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	float oneOverd = tan(m_camera.FOV * PI / (2 * 180));
	GzMatrix Xpi_temp =
	{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, oneOverd, 0,
		0, 0, oneOverd, 1
	};

	memmove(m_camera.Xpi, Xpi_temp, sizeof(m_camera.Xpi));

	//Xiw
	GzCoord cl =
	{
		m_camera.lookat[0] - m_camera.position[0],
		m_camera.lookat[1] - m_camera.position[1],
		m_camera.lookat[2] - m_camera.position[2],
	};
	GzCoord z_result =
	{
		cl[0] / calculateNorm(cl),
		cl[1] / calculateNorm(cl),
		cl[2] / calculateNorm(cl)
	};
	float temp = calculateDotProduct(m_camera.worldup, z_result);
	GzCoord toSub =
	{
		temp*z_result[0],
		temp*z_result[1],
		temp*z_result[2]
	};
	GzCoord upPrime =
	{
		m_camera.worldup[0] - toSub[0],
		m_camera.worldup[1] - toSub[1],
		m_camera.worldup[2] - toSub[2]
	};
	GzCoord y_result =
	{
		upPrime[0] / calculateNorm(upPrime),
		upPrime[1] / calculateNorm(upPrime),
		upPrime[2] / calculateNorm(upPrime)
	};
	GzCoord x_result;
	calculateCrossProduct(y_result, z_result, &x_result[0], &x_result[1], &x_result[2]);

	GzMatrix Xiw_temp =
	{
		x_result[0], x_result[1], x_result[2], -calculateDotProduct(x_result, m_camera.position),
		y_result[0], y_result[1], y_result[2], -calculateDotProduct(y_result, m_camera.position),
		z_result[0], z_result[1], z_result[2], -calculateDotProduct(z_result, m_camera.position),
		0, 0, 0, 1
	};
	memmove(m_camera.Xiw, Xiw_temp, sizeof(m_camera.Xiw));

	GzPushMatrix(Xsp);

	GzPushMatrix(m_camera.Xpi);

	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/
	m_camera = camera;
	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	GzMatrix I =
	{
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1
	};

	if (matlevel == 0)
	{
		memmove(Ximage[matlevel], matrix, sizeof(Ximage[matlevel]));
		GzPushMatrixOnNormStack(I);
		matlevel = matlevel + 1;
		return GZ_SUCCESS;
	}
	else if (matlevel == 1)
	{
		GzMatrix productMatrix;
		multiplyMatrices(Ximage[matlevel - 1], matrix, &productMatrix);
		memmove(Ximage[matlevel], productMatrix, sizeof(Ximage[matlevel]));
		GzPushMatrixOnNormStack(I);
		matlevel = matlevel + 1;
		return GZ_SUCCESS;
	}
	else if (matlevel == MATLEVELS)
	{
		return GZ_FAILURE;		//stack full
	}
	else
	{
		GzMatrix productMatrix;
		multiplyMatrices(Ximage[matlevel - 1], matrix, &productMatrix);
		memmove(Ximage[matlevel], productMatrix, sizeof(Ximage[matlevel]));
		GzPushMatrixOnNormStack(matrix);
		matlevel = matlevel + 1;
		return GZ_SUCCESS;
	}
}

int GzRender::GzPushMatrixOnNormStack(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	if ((matlevel) == MATLEVELS)
	{
		return GZ_FAILURE;		//stack full
	}

	// get rid of translation
	matrix[0][3] = 0;
	matrix[1][3] = 0;
	matrix[2][3] = 0;

	GzCoord firstRow = { matrix[0][0], matrix[0][1], matrix[0][2] };
	float k = calculateNorm(firstRow);

	//normalise unitary rotation
	if (k != 0)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				matrix[i][j] = matrix[i][j] / k;
			}
		}
	}


	if ((matlevel) == 0)
	{
		memmove(Xnorm[(matlevel)], matrix, sizeof(Xnorm[(matlevel)]));
		return GZ_SUCCESS;
	}
	else
	{
		GzMatrix productMatrix;
		multiplyMatrices(Xnorm[matlevel - 1], matrix, &productMatrix);
		memmove(Xnorm[matlevel], productMatrix, sizeof(Xnorm[matlevel]));
		return GZ_SUCCESS;
	}
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel <= 0)
		return GZ_FAILURE;		// empty stack
	matlevel = matlevel - 1;

	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */

	// do nothing if pixel is not on screen
	if (i >= 0 && i < xres && j >= 0 && j < yres)
	{
		//clamp color values if out of range
		GzPixel current;
		if (r < 0)
			current.red = 0;
		if (r > 4095)
			current.red = 4095;
		else
			current.red = r;

		if (g < 0)
			current.green = 0;
		if (g > 4095)
			current.green = 4095;
		else
			current.green = g;

		if (b < 0)
			current.blue = 0;
		if (b > 4095)
			current.blue = 4095;
		else
			current.blue = b;

		if (a < 0)
			current.alpha = 0;
		if (a > 4095)
			current.alpha = 4095;
		else
			current.alpha = a;

		current.z = z;
		pixelbuffer[j*xres + i] = current;		// add the pixel to correct location in array
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i >= 0 && i < xres && j >= 0 && j < yres)
	{
		GzPixel fetched = pixelbuffer[j*xres + i];
		*r = fetched.red;
		*g = fetched.green;
		*b = fetched.blue;
		*a = fetched.alpha;
		*z = fetched.z;
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	GzIntensity r, g, b, a;
	int rWrite, gWrite, bWrite;
	GzDepth z;
	fprintf(outfile, "P6 %d %d 255\r", xres, yres);		// did not work with P6
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			GzGet(i, j, &r, &g, &b, &a, &z);
			// short to 8 bit conversion
			rWrite = r >> 4;
			gWrite = g >> 4;
			bWrite = b >> 4;
			fprintf(outfile, "%c%c%c", rWrite, gWrite, bWrite);
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
	- NOT red, green, and blue !!!
	*/
	GzIntensity r, g, b, a;
	GzDepth z;
	int index = 0;
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			GzGet(i, j, &r, &g, &b, &a, &z);
			framebuffer[index] = (b >> 4);
			framebuffer[index + 1] = (g >> 4);
			framebuffer[index + 2] = (r >> 4);
			index = index + 3;
		}
	}
	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	/*
	- GzPutAttribute() must accept the following tokens/values:

	- GZ_RGB_COLOR					GzColor		default flat-shade color
	- GZ_INTERPOLATE				int			shader interpolation mode
	- GZ_DIRECTIONAL_LIGHT			GzLight
	- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
	- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
	- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
	- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
	- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
	*/
	for (int i = 0; i < numAttributes; i++)
	{
		//if token name is GZ_RGB_COLOR put token value in flatcolor
		switch (nameList[i])
		{
		case GZ_RGB_COLOR:
			for (int j = 0; j < 3; j++)
			{
				flatcolor[j] = (*(GzColor*)valueList[i])[j];
			}
			break;

		case GZ_INTERPOLATE:
			interp_mode = *(int *)valueList[i];
			break;

		case GZ_DIRECTIONAL_LIGHT:
			lights[numlights] = *(GzLight*)valueList[i];
			numlights = numlights + 1;
			break;

		case GZ_AMBIENT_LIGHT:
			ambientlight = *(GzLight*)valueList[i];
			break;

		case GZ_DIFFUSE_COEFFICIENT:
			for (int j = 0; j < 3; j++)
			{
				Kd[j] = (*(GzColor*)valueList[i])[j];
			}
			break;

		case GZ_AMBIENT_COEFFICIENT:
			for (int j = 0; j < 3; j++)
			{
				Ka[j] = (*(GzColor*)valueList[i])[j];
			}
			break;

		case GZ_SPECULAR_COEFFICIENT:
			for (int j = 0; j < 3; j++)
			{
				Ks[j] = (*(GzColor*)valueList[i])[j];
			}
			break;

		case GZ_DISTRIBUTION_COEFFICIENT:
			spec = *(float*)valueList[i];
			break;

		case GZ_TEXTURE_MAP:
			tex_fun = (GzTexture)valueList[i];
			break;
		case GZ_REND_NO:
			rendIndex = *(int*)valueList[i];
			break;
		case GZ_AASHIFTX:
			offsetX[rendIndex] = *(float*)valueList[i];
			break;
		case GZ_AASHIFTY:
			offsetY[rendIndex] = *(float*)valueList[i];
			break;
		}
	}
	return GZ_SUCCESS;
}

// Helper methods needed for putTriangle
bool GzRender::isHorizontal(float y0, float y1)
{
	if (y0 == y1)
		return TRUE;
	return FALSE;
}

void GzRender::computeBoundingBox(float x0, float x1, float x2, float y0, float y1, float y2,
	float *xmin, float *xmax, float *ymin, float *ymax)
{
	if (x0 < x1) {
		if (x2 < x0)
			*xmin = x2;
		else
			*xmin = x0;
	}
	else {
		if (x2 < x1)
			*xmin = x2;
		else
			*xmin = x1;
	}

	if (x0 > x1) {
		if (x2 > x0)
			*xmax = x2;
		else
			*xmax = x0;
	}
	else {
		if (x2 > x1)
			*xmax = x2;
		else
			*xmax = x1;
	}

	if (y0 < y1) {
		if (y2 < y0)
			*ymin = y2;
		else
			*ymin = y0;
	}
	else {
		if (y2 < y1)
			*ymin = y2;
		else
			*ymin = y1;
	}

	if (y0 > y1) {
		if (y2 > y0)
			*ymax = y2;
		else
			*ymax = y0;
	}
	else {
		if (y2 > y1)
			*ymax = y2;
		else
			*ymax = y1;
	}
}

void GzRender::computeEquation(float tailX, float tailY, float headX, float headY,
	float *a, float *b, float *c)
{
	float dy = headY - tailY;
	float dx = headX - tailX;
	*a = dy;
	*b = -dx;
	*c = (tailY * dx) - (tailX * dy);
}

int GzRender::evaluateEquation(float a, float b, float c, float x, float y)
{
	float val = (a*x + b*y + c);
	if (val > 0)
		return 1;
	else if (val < 0)
		return -1;
	else
		return 0;
}

int GzRender::assignFlavorsHorizontal(GzCoord vertices[], int hor_v1, int hor_v2, int other,
	int *e01_flavor, int *e12_flavor, int *e20_flavor)
{
	// e01 is horizontal
	if (hor_v1 == 0 && hor_v2 == 1)
	{
		// compare y values of horizontal edge and other
		if (vertices[other][1] < vertices[hor_v1][1])
		{
			*e01_flavor = TOP;
		}
		else
		{
			*e01_flavor = BOTTOM;
		}

		// the left most one with other forms L edge
		if (vertices[0][0] < vertices[1][0])
		{
			*e20_flavor = LEFT;
			*e12_flavor = RIGHT;
		}
		else
		{
			*e20_flavor = RIGHT;
			*e12_flavor = LEFT;
		}
		return GZ_SUCCESS;
	}

	// e12 is horizontal
	if (hor_v1 == 1 && hor_v2 == 2)
	{
		// compare y values of horizontal edge and other
		if (vertices[other][1] < vertices[hor_v1][1])
		{
			*e12_flavor = TOP;
		}
		else
		{
			*e12_flavor = BOTTOM;
		}

		// the left most one with other forms L edge
		if (vertices[1][0] < vertices[2][0])
		{
			*e01_flavor = LEFT;
			*e20_flavor = RIGHT;
		}
		else
		{
			*e01_flavor = RIGHT;
			*e20_flavor = LEFT;
		}
		return GZ_SUCCESS;
	}

	// e20 is horizontal
	if (hor_v1 == 2 && hor_v2 == 0)
	{
		// compare y values of horizontal edge and other
		if (vertices[other][1] < vertices[hor_v1][1])
		{
			*e20_flavor = TOP;
		}
		else
		{
			*e20_flavor = BOTTOM;
		}

		// the left most one with other forms L edge
		if (vertices[2][0] < vertices[0][0])
		{
			*e12_flavor = LEFT;
			*e01_flavor = RIGHT;
		}
		else
		{
			*e12_flavor = RIGHT;
			*e01_flavor = LEFT;
		}
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}

int GzRender::assignFlavorsGeneral(GzCoord vertices[], float ymin, float ymax,
	int *e01_flavor, int *e12_flavor, int *e20_flavor)
{
	int other, maxV, minV;			// indices of vertices which are max, min and other
	for (int i = 0; i < 3; i++)
	{
		if (vertices[i][1] == ymin)
		{
			minV = i;
		}
		else if (vertices[i][1] == ymax)
		{
			maxV = i;
		}
		else
		{
			other = i;
		}
	}

	float a, b, c;			// co-efficients of equation of min-max edge
	computeEquation(vertices[minV][0], vertices[minV][1], vertices[maxV][0], vertices[maxV][1],
		&a, &b, &c);

	float x_comp = vertices[other][0];
	float y_comp = vertices[other][1];
	float x_calc = ((b*y_comp) + c) / (-a);

	// e01 is sole edge
	if ((minV == 0 && maxV == 1) || (minV == 1 && maxV == 0))
	{
		if (x_calc < x_comp)
		{
			*e01_flavor = LEFT;
			*e12_flavor = RIGHT;
			*e20_flavor = RIGHT;
			return GZ_SUCCESS;
		}
		else
		{
			*e01_flavor = RIGHT;
			*e12_flavor = LEFT;
			*e20_flavor = LEFT;
			return GZ_SUCCESS;
		}
	}

	// e12 is sole edge
	if ((minV == 2 && maxV == 1) || (minV == 1 && maxV == 2))
	{
		if (x_calc < x_comp)
		{
			*e12_flavor = LEFT;
			*e01_flavor = RIGHT;
			*e20_flavor = RIGHT;
			return GZ_SUCCESS;
		}
		else
		{
			*e12_flavor = RIGHT;
			*e01_flavor = LEFT;
			*e20_flavor = LEFT;
			return GZ_SUCCESS;
		}
	}

	// e20 is sole edge
	if ((minV == 0 && maxV == 2) || (minV == 2 && maxV == 0))
	{
		if (x_calc < x_comp)
		{
			*e20_flavor = LEFT;
			*e12_flavor = RIGHT;
			*e01_flavor = RIGHT;
			return GZ_SUCCESS;
		}
		else
		{
			*e20_flavor = RIGHT;
			*e12_flavor = LEFT;
			*e01_flavor = LEFT;
			return GZ_SUCCESS;
		}
	}

	return GZ_FAILURE;
}

void GzRender::getCWEdgesHori(GzCoord vertices[], int hori_v1, int hori_v2, int other, int hori_flavor,
	int *e1_head, int *e1_tail, int *e2_head, int *e2_tail, int *e3_head, int *e3_tail)
{
	if (hori_flavor == TOP)
	{
		if (vertices[hori_v1][0] < vertices[hori_v2][0])
		{
			*e1_tail = hori_v1;
			*e1_head = hori_v2;
			*e2_tail = hori_v2;
			*e2_head = other;
			*e3_tail = other;
			*e3_head = hori_v1;
		}
		else
		{
			*e1_tail = hori_v2;
			*e1_head = hori_v1;
			*e2_tail = hori_v1;
			*e2_head = other;
			*e3_tail = other;
			*e3_head = hori_v2;
		}
	}
	else
		// BOTTOM case
	{
		if (vertices[hori_v1][0] > vertices[hori_v2][0])
		{
			*e1_tail = hori_v1;
			*e1_head = hori_v2;
			*e2_tail = hori_v2;
			*e2_head = other;
			*e3_tail = other;
			*e3_head = hori_v1;
		}
		else
		{
			*e1_tail = hori_v2;
			*e1_head = hori_v1;
			*e2_tail = hori_v1;
			*e2_head = other;
			*e3_tail = other;
			*e3_head = hori_v2;
		}
	}
}

void GzRender::getCWEdgesGeneral(GzCoord vertices[], int soleEdge_flavor,
	int minV, int maxV, int midV, int *e1_head, int *e1_tail, int *e1_flavor,
	int *e2_head, int *e2_tail, int *e2_flavor, int *e3_head, int *e3_tail, int *e3_flavor)
{
	if (soleEdge_flavor == RIGHT)
	{
		*e1_tail = minV;
		*e1_head = maxV;
		*e1_flavor = RIGHT;
		*e2_tail = maxV;
		*e2_head = midV;
		*e2_flavor = LEFT;
		*e3_tail = midV;
		*e3_head = minV;
		*e3_flavor = LEFT;
	}
	else if (soleEdge_flavor == LEFT)
	{
		*e1_tail = maxV;
		*e1_head = minV;
		*e1_flavor = LEFT;
		*e2_tail = minV;
		*e2_head = midV;
		*e2_flavor = RIGHT;
		*e3_tail = midV;
		*e3_head = maxV;
		*e3_flavor = RIGHT;
	}
}

//4x4 from stack * 4x1 vertex
void GzRender::applyTransformation(GzMatrix stack, float vertex[4],
	float *resultx, float *resulty, float *resultz, float *resultw)
{
	float homogeneousResult[4];
	for (int i = 0; i < 4; i++)
	{
		float sum = 0;
		for (int j = 0; j < 4; j++)
		{
			sum = sum + (stack[i][j] * vertex[j]);
		}
		homogeneousResult[i] = sum;
	}
	*resultx = homogeneousResult[0];
	*resulty = homogeneousResult[1];
	*resultz = homogeneousResult[2];
	*resultw = homogeneousResult[3];
}

int GzRender::checkNormalSide(GzCoord N, GzCoord L, GzCoord E, float *nx, float *ny, float *nz)
{
	*nx = N[0];
	*ny = N[1];
	*nz = N[2];

	if (calculateDotProduct(N, L) < 0 && calculateDotProduct(N, E) < 0)
	{
		*nx = -1 * N[0];
		*ny = -1 * N[1];
		*nz = -1 * N[2];
		return -1;
	}
	if (calculateDotProduct(N, L) > 0 && calculateDotProduct(N, E) > 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}

}

void GzRender::computeShadingEquation(GzCoord transformedNormal, float *r, float *g, float *b,
	GzColor Ka_tex, GzColor Kd_tex, GzColor Ks_tex)
{
	GzColor computedColor;
	GzCoord E = { 0, 0, -1 };

	GzColor ambientTerm, specTerm, diffTerm;

	//initialise all to 0 to avoid garbage values in sum
	for (int j = 0; j < 3; j++)
	{
		ambientTerm[j] = 0;
		specTerm[j] = 0;
		diffTerm[j] = 0;
	}

	//each color
	for (int j = 0; j < 3; j++)
	{
		ambientTerm[j] = Ka_tex[j] * ambientlight.color[j];
		//each light
		for (int k = 0; k < numlights; k++)
		{
			float nx, ny, nz;
			int result = checkNormalSide(transformedNormal, lights[k].direction, E, &nx, &ny, &nz);
			if (result == 0)
				continue;
			GzCoord finalNormal;
			finalNormal[0] = nx;
			finalNormal[1] = ny;
			finalNormal[2] = nz;
			float dotProduct = calculateDotProduct(finalNormal, lights[k].direction);
			diffTerm[j] = diffTerm[j] + (lights[k].color[j] * dotProduct);

			GzCoord R;
			float tempX, tempY, tempZ;
			//2(N.L)N
			scalarVectorMultiply(2 * dotProduct, finalNormal, &tempX, &tempY, &tempZ);
			R[0] = tempX;
			R[1] = tempY;
			R[2] = tempZ;
			//R = 2(N.L)N - L
			calculateDifferenceVectors(R, lights[k].direction, &tempX, &tempY, &tempZ);
			R[0] = tempX;
			R[1] = tempY;
			R[2] = tempZ;

			//normalise
			float length = calculateNorm(R);
			R[0] = R[0] / length;
			R[1] = R[1] / length;
			R[2] = R[2] / length;

			dotProduct = calculateDotProduct(R, E);
			//clamp to avoid negative R.E case
			if (dotProduct < 0)
				dotProduct = 0;
			specTerm[j] = specTerm[j] + (lights[k].color[j] * pow(dotProduct, spec));
		}

		diffTerm[j] = Kd_tex[j] * diffTerm[j];
		specTerm[j] = Ks_tex[j] * specTerm[j];

		computedColor[j] = specTerm[j] + diffTerm[j] + ambientTerm[j];
		//clamping to avoid overflow
		if (computedColor[j] > 1)
			computedColor[j] = 1;
		if (computedColor[j] < 0)
			computedColor[j] = 0;
	}

	*r = computedColor[0];
	*g = computedColor[1];
	*b = computedColor[2];
}

float GzRender::interpolate(GzCoord v1, GzCoord v2, GzCoord v3, int x, int y)
{
	float tempx, tempy, tempz;
	GzCoord e1, e2;
	float a, b, c, d;
	float result;

	calculateDifferenceVectors(v1, v2, &tempx, &tempy, &tempz);
	e1[0] = tempx;
	e1[1] = tempy;
	e1[2] = tempz;

	calculateDifferenceVectors(v3, v2, &tempx, &tempy, &tempz);
	e2[0] = tempx;
	e2[1] = tempy;
	e2[2] = tempz;

	calculateCrossProduct(e1, e2, &a, &b, &c);
	d = -((a * v1[0]) + (b * v1[1]) + (c * v1[2]));

	result = -((a * x) + (b * y) + d) / c;

	return result;
}

void GzRender::rasterize(GzCoord transformedVertices[], GzCoord transformedNormals[], GzTextureIndex texels[])
{
	float x1, x2, x0, y1, y2, y0, z1, z2, z0;
	x0 = transformedVertices[0][0];
	y0 = transformedVertices[0][1];
	z0 = transformedVertices[0][2];

	x1 = transformedVertices[1][0];
	y1 = transformedVertices[1][1];
	z1 = transformedVertices[1][2];

	x2 = transformedVertices[2][0];
	y2 = transformedVertices[2][1];
	z2 = transformedVertices[2][2];

	/*
	1. Check if any horizontal edge
	1.1 Yes? Assign T/B, L, R to all 3 edges
	1.2 No? Assign L, R to all 3 edges

	2. Assign CW edges (head/tail)
	2.1 Horizontal edge
	2.1.1 T? E1: Xmin -> Xmax, E2: Xmax -> other, E3: other -> Xmin
	2.1.2 B? E1: Xmax -> Xmin, E2: Xmin -> other, E3: other -> Xmax
	2.2 No Horizontal edge, check flavor of Vmin-Vmax edge
	2.2.1 R? E1: Vmin -> Vmax, E2: Vmax -> other, E3: other -> Vmin
	2.2.2 L? E1: Vmax -> Vmin, E2: Vmin -> other, E3: other -> Vmax

	3. Compute bounding box: For full tri, get Xmin, Xmax, Ymin, Ymax

	4. Compute Z interpolation A, B, C, D

	5. Scale flatcolor - get r, g, b

	5. Double loop through BB
	5.1 Check if Ax + By + C > 0 for all 3 edges
	5.2 Check if Ax + By + C = 0 for any one edge
	5.3 Interpolate Z for that pixel
	5.4 Call Gzput with (x, y, r, g, b, 1, z_inter)
	*/

	bool hasHorizontalEdge = FALSE;
	int hori_v1 = 0;
	int hori_v2 = 1;
	int other = 2;
	if (isHorizontal(y0, y1) || isHorizontal(y1, y2) || isHorizontal(y0, y2))
	{
		hasHorizontalEdge = TRUE;
		if (isHorizontal(y0, y1))
		{
			hori_v1 = 0;
			hori_v2 = 1;
			other = 2;
		}
		else if (isHorizontal(y1, y2))
		{
			hori_v1 = 1;
			hori_v2 = 2;
			other = 0;
		}
		else if (isHorizontal(y0, y2))
		{
			hori_v1 = 2;
			hori_v2 = 0;
			other = 1;
		}
	}

	int e01_flavor, e12_flavor, e20_flavor;
	float xmin, xmax, ymin, ymax;

	computeBoundingBox(x0, x1, x2, y0, y1, y2, &xmin, &xmax, &ymin, &ymax);

	if (hasHorizontalEdge)
		assignFlavorsHorizontal(transformedVertices, hori_v1, hori_v2, other,
			&e01_flavor, &e12_flavor, &e20_flavor);
	else
		assignFlavorsGeneral(transformedVertices, ymin, ymax,
			&e01_flavor, &e12_flavor, &e20_flavor);

	int e1_head, e1_tail, e1_flavor;
	int e2_head, e2_tail, e2_flavor;
	int e3_head, e3_tail, e3_flavor;

	if (hasHorizontalEdge)
	{
		if (hori_v1 == 0 && hori_v2 == 1)
		{
			getCWEdgesHori(transformedVertices, hori_v1, hori_v2, other, e01_flavor,
				&e1_head, &e1_tail, &e2_head, &e2_tail, &e3_head, &e3_tail);
		}
		else if (hori_v1 == 1 && hori_v2 == 2)
		{
			getCWEdgesHori(transformedVertices, hori_v1, hori_v2, other, e12_flavor,
				&e1_head, &e1_tail, &e2_head, &e2_tail, &e3_head, &e3_tail);
		}
		else if (hori_v1 == 2 && hori_v2 == 0)
		{
			getCWEdgesHori(transformedVertices, hori_v1, hori_v2, other, e20_flavor,
				&e1_head, &e1_tail, &e2_head, &e2_tail, &e3_head, &e3_tail);
		}
	}
	else
	{
		int midV, maxV, minV;			// indices of transformedVertices which are max, min and other
		for (int i = 0; i < 3; i++)
		{
			if (transformedVertices[i][1] == ymin)
			{
				minV = i;
			}
			else if (transformedVertices[i][1] == ymax)
			{
				maxV = i;
			}
			else
			{
				midV = i;
			}
		}

		int soleEdge_flavor;

		if ((minV == 0 && maxV == 1) || (minV == 1 && maxV == 0))
		{
			soleEdge_flavor = e01_flavor;
		}
		else if ((minV == 2 && maxV == 1) || (minV == 1 && maxV == 2))
		{
			soleEdge_flavor = e12_flavor;
		}
		else if ((minV == 0 && maxV == 2) || (minV == 2 && maxV == 0))
		{
			soleEdge_flavor = e20_flavor;
		}
		getCWEdgesGeneral(transformedVertices, soleEdge_flavor,
			minV, maxV, midV, &e1_head, &e1_tail, &e1_flavor,
			&e2_head, &e2_tail, &e2_flavor, &e3_head, &e3_tail, &e3_flavor);
	}

	// Z interpolation - compute A, B, C, D
	float a_z, b_z, c_z, d_z;
	float e0_x, e1_x, e0_y, e1_y, e0_z, e1_z;		// edge vectors
	e0_x = transformedVertices[e1_head][0] - transformedVertices[e1_tail][0];
	e0_y = transformedVertices[e1_head][1] - transformedVertices[e1_tail][1];
	e0_z = transformedVertices[e1_head][2] - transformedVertices[e1_tail][2];

	e1_x = transformedVertices[e2_head][0] - transformedVertices[e2_tail][0];
	e1_y = transformedVertices[e2_head][1] - transformedVertices[e2_tail][1];
	e1_z = transformedVertices[e2_head][2] - transformedVertices[e2_tail][2];

	a_z = (e0_y * e1_z) - (e0_z * e1_y);
	b_z = (e0_z * e1_x) - (e0_x * e1_z);
	c_z = (e0_x * e1_y) - (e0_y * e1_x);

	d_z = -((a_z * transformedVertices[0][0]) + (b_z * transformedVertices[0][1]) + (c_z * transformedVertices[0][2]));

	// co-effs for the edge equations
	float e1_a, e1_b, e1_c;
	float e2_a, e2_b, e2_c;
	float e3_a, e3_b, e3_c;

	computeEquation(transformedVertices[e1_tail][0], transformedVertices[e1_tail][1],
		transformedVertices[e1_head][0], transformedVertices[e1_head][1], &e1_a, &e1_b, &e1_c);
	computeEquation(transformedVertices[e2_tail][0], transformedVertices[e2_tail][1],
		transformedVertices[e2_head][0], transformedVertices[e2_head][1], &e2_a, &e2_b, &e2_c);
	computeEquation(transformedVertices[e3_tail][0], transformedVertices[e3_tail][1],
		transformedVertices[e3_head][0], transformedVertices[e3_head][1], &e3_a, &e3_b, &e3_c);

	// perspective texel transformation
	GzTextureIndex perspectiveTexels[3];
	for (int i = 0; i < 3; i++)
	{
		float screenZ = transformedVertices[i][2];
		float factor = screenZ / ((float)MAXINT - screenZ);
		perspectiveTexels[i][0] = texels[i][0] / (factor + 1);
		perspectiveTexels[i][1] = texels[i][1] / (factor + 1);
	}

	float z;
	GzColor vertexColors[3];

	// do once per triangle
	if (interp_mode == GZ_COLOR)
	{
		//no texture
		if (tex_fun == 0)
		{
			//for each vertex
			for (int i = 0; i < 3; i++)
			{
				float r, g, b;
				computeShadingEquation(transformedNormals[i], &r, &g, &b, Ka, Kd, Ks);
				vertexColors[i][0] = r;
				vertexColors[i][1] = g;
				vertexColors[i][2] = b;
			}
		}
		else
		{
			//for each vertex
			for (int i = 0; i < 3; i++)
			{
				float r, g, b;
				GzColor Ka_temp = { 1, 1, 1 };
				GzColor Kd_temp = { 1, 1, 1 };
				GzColor Ks_temp = { 1, 1, 1 };
				computeShadingEquation(transformedNormals[i], &r, &g, &b, Ka_temp, Kd_temp, Ks_temp);
				vertexColors[i][0] = r;
				vertexColors[i][1] = g;
				vertexColors[i][2] = b;
			}
		}
	}

	for (int x = xmin; x <= xmax; x++)
	{
		for (int y = ymin; y <= ymax; y++)
		{
			int result1 = evaluateEquation(e1_a, e1_b, e1_c, x, y);
			int result2 = evaluateEquation(e2_a, e2_b, e2_c, x, y);
			int result3 = evaluateEquation(e3_a, e3_b, e3_c, x, y);

			if (
				((result1 == 0) && ((e1_flavor == TOP) || (e1_flavor == LEFT))) ||
				((result2 == 0) && ((e2_flavor == TOP) || (e2_flavor == LEFT))) ||
				((result3 == 0) && ((e3_flavor == TOP) || (e3_flavor == LEFT))) ||
				((result1 < 0) && (result2 < 0) && (result3 < 0))
				)
			{
				z = -((a_z * x) + (b_z * y) + d_z) / c_z;
				GzDepth z_fetch;
				GzIntensity r_fetch, g_fetch, b_fetch, a_fetch;
				GzGet(x, y, &r_fetch, &g_fetch, &b_fetch, &a_fetch, &z_fetch);

				if (z < z_fetch)
				{
					GzTextureIndex transformedUV;

					// interpolation for u & v in perspective space
					for (int i = 0; i < 2; i++)
					{
						GzCoord tex1, tex2, tex3;
						//tex1[0] = transformedVertices[0][0];		//x
						//tex1[1] = transformedVertices[0][1];		//y
						//tex1[2] = texels[0][i];			//u

						//tex2[0] = transformedVertices[1][0];		//x
						//tex2[1] = transformedVertices[1][1];		//y
						//tex2[2] = texels[1][i];			//u

						//tex3[0] = transformedVertices[2][0];		//x
						//tex3[1] = transformedVertices[2][1];		//y
						//tex3[2] = texels[2][i];			//u


						tex1[0] = transformedVertices[0][0];		//x
						tex1[1] = transformedVertices[0][1];		//y
						tex1[2] = perspectiveTexels[0][i];			//u,v

						tex2[0] = transformedVertices[1][0];		//x
						tex2[1] = transformedVertices[1][1];		//y
						tex2[2] = perspectiveTexels[1][i];			//u,v

						tex3[0] = transformedVertices[2][0];		//x
						tex3[1] = transformedVertices[2][1];		//y
						tex3[2] = perspectiveTexels[2][i];			//u,v

						transformedUV[i] = interpolate(tex1, tex2, tex3, x, y);
					}

					// get back to affine space
					GzTextureIndex uv;

					float back_factor = z / ((float)MAXINT - z);
					uv[0] = transformedUV[0] * (back_factor + 1);
					uv[1] = transformedUV[1] * (back_factor + 1);

					//flat color
					if (interp_mode == GZ_FLAT)
					{
						float r, g, b;
						computeShadingEquation(transformedNormals[0], &r, &g, &b, Ka, Kd, Ks);		//color of first vertex
						GzPut(x, y, ctoi(r), ctoi(g), ctoi(b), 1, z);
					}
					//Gouraud shading
					else if (interp_mode == GZ_COLOR)
					{
						//interpolate R,G,B from vertex colors and call GzPut

						//interpolating colors
						GzColor finalColor;
						for (int i = 0; i < 3; i++)
						{
							GzCoord v1, v2, v3;
							v1[0] = transformedVertices[0][0];		//x
							v1[1] = transformedVertices[0][1];		//y
							v1[2] = vertexColors[0][i];				//r,g,b

							v2[0] = transformedVertices[1][0];		//x
							v2[1] = transformedVertices[1][1];		//y
							v2[2] = vertexColors[1][i];				//r,g,b

							v3[0] = transformedVertices[2][0];		//x
							v3[1] = transformedVertices[2][1];		//y
							v3[2] = vertexColors[2][i];				//r,g,b

							if (tex_fun == 0)
							{
								finalColor[i] = interpolate(v1, v2, v3, x, y);
							}
							else
							{
								GzColor texColor;
								tex_fun(uv[0], uv[1], texColor);
								//tex_fun(transformedUV[0], transformedUV[1], texColor);
								finalColor[i] = texColor[i] * interpolate(v1, v2, v3, x, y);
							}
						}

						GzPut(x, y, ctoi(finalColor[0]), ctoi(finalColor[1]), ctoi(finalColor[2]), 1, z);
					}
					//Phong shading
					else if (interp_mode == GZ_NORMALS)
					{
						// interpolate Nx, Ny, Nz from transformedNormals
						// compute shading equation using that normal
						// call gzput using r, g, b returned
						GzCoord interpolatedNormal;
						//x, y, z
						for (int i = 0; i < 3; i++)
						{
							GzCoord v1, v2, v3;
							v1[0] = transformedVertices[0][0];		//x
							v1[1] = transformedVertices[0][1];		//y
							v1[2] = transformedNormals[0][i];		//n component

							v2[0] = transformedVertices[1][0];		//x
							v2[1] = transformedVertices[1][1];		//y
							v2[2] = transformedNormals[1][i];		//n component

							v3[0] = transformedVertices[2][0];		//x
							v3[1] = transformedVertices[2][1];		//y
							v3[2] = transformedNormals[2][i];		//n component

							interpolatedNormal[i] = interpolate(v1, v2, v3, x, y);
						}

						//normalise
						float length = calculateNorm(interpolatedNormal);
						interpolatedNormal[0] = interpolatedNormal[0] / length;
						interpolatedNormal[1] = interpolatedNormal[1] / length;
						interpolatedNormal[2] = interpolatedNormal[2] / length;

						float r, g, b;
						if (tex_fun == 0)
						{
							computeShadingEquation(interpolatedNormal, &r, &g, &b, Ka, Kd, Ks);
						}
						else
						{
							GzColor Kt;
							tex_fun(uv[0], uv[1], Kt);
							//tex_fun(transformedUV[0], transformedUV[1], texColor);
							computeShadingEquation(interpolatedNormal, &r, &g, &b, Kt, Kt, Ks);
						}

						GzPut(x, y, ctoi(r), ctoi(g), ctoi(b), 1, z);
					}
				}
			}
		}
	}
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList, int rendNo)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
	GZ_NULL_TOKEN:		do nothing - no values
	GZ_POSITION:		3 vert positions in model space
	-- Return error code
	*/
	/*
	-- Xform positions of verts using matrix on top of stack
	-- Clip - just discard any triangle with any vert(s) behind view plane
	- optional: test for triangles with all three verts off-screen (trivial frustum cull)
	-- invoke triangle rasterizer
	*/
	GzCoord vertices[3];
	GzCoord normals[3];
	GzTextureIndex texels[3];		// (u,v) per vertex
	for (int i = 0; i < numParts; i++)
	{
		if (nameList[i] == GZ_POSITION)
		{
			GzCoord *vertexPointer = (GzCoord*)valueList[i];

			vertices[0][0] = vertexPointer[0][0];
			vertices[0][1] = vertexPointer[0][1];
			vertices[0][2] = vertexPointer[0][2];

			vertices[1][0] = vertexPointer[1][0];
			vertices[1][1] = vertexPointer[1][1];
			vertices[1][2] = vertexPointer[1][2];

			vertices[2][0] = vertexPointer[2][0];
			vertices[2][1] = vertexPointer[2][1];
			vertices[2][2] = vertexPointer[2][2];
		}

		if (nameList[i] == GZ_NORMAL)
		{
			GzCoord *normalPointer = (GzCoord*)valueList[i];

			normals[0][0] = normalPointer[0][0];
			normals[0][1] = normalPointer[0][1];
			normals[0][2] = normalPointer[0][2];

			normals[1][0] = normalPointer[1][0];
			normals[1][1] = normalPointer[1][1];
			normals[1][2] = normalPointer[1][2];

			normals[2][0] = normalPointer[2][0];
			normals[2][1] = normalPointer[2][1];
			normals[2][2] = normalPointer[2][2];
		}

		if (nameList[i] == GZ_TEXTURE_INDEX)
		{
			GzTextureIndex *uvPointer = (GzTextureIndex*)valueList[i];

			texels[0][0] = uvPointer[0][0];
			texels[0][1] = uvPointer[0][1];

			texels[1][0] = uvPointer[1][0];
			texels[1][1] = uvPointer[1][1];

			texels[2][0] = uvPointer[2][0];
			texels[2][1] = uvPointer[2][1];
		}

		if (nameList[i] == GZ_NULL_TOKEN)
		{
			return GZ_FAILURE;
		}
	}

	float homogenousV0[4] =
	{
		vertices[0][0],
		vertices[0][1],
		vertices[0][2],
		1
	};
	float homogenousV1[4] =
	{
		vertices[1][0],
		vertices[1][1],
		vertices[1][2],
		1
	};
	float homogenousV2[4] =
	{
		vertices[2][0],
		vertices[2][1],
		vertices[2][2],
		1
	};

	GzCoord transformedVertices[3];
	float resultX, resultY, resultZ, resultW;

	applyTransformation(Ximage[matlevel - 1], homogenousV0, &resultX, &resultY, &resultZ, &resultW);
	if (resultZ < 0)
		return GZ_SUCCESS;

	transformedVertices[0][0] = resultX / resultW;
	transformedVertices[0][1] = resultY / resultW;
	transformedVertices[0][2] = resultZ / resultW;

	applyTransformation(Ximage[matlevel - 1], homogenousV1, &resultX, &resultY, &resultZ, &resultW);
	if (resultZ < 0)
		return GZ_SUCCESS;

	transformedVertices[1][0] = resultX / resultW;
	transformedVertices[1][1] = resultY / resultW;
	transformedVertices[1][2] = resultZ / resultW;

	applyTransformation(Ximage[matlevel - 1], homogenousV2, &resultX, &resultY, &resultZ, &resultW);
	if (resultZ < 0)
		return GZ_SUCCESS;

	transformedVertices[2][0] = resultX / resultW;
	transformedVertices[2][1] = resultY / resultW;
	transformedVertices[2][2] = resultZ / resultW;

	GzCoord transformedNormals[3];
	float homogenousN0[4] =
	{
		normals[0][0],
		normals[0][1],
		normals[0][2],
		1
	};
	float homogenousN1[4] =
	{
		normals[1][0],
		normals[1][1],
		normals[1][2],
		1
	};
	float homogenousN2[4] =
	{
		normals[2][0],
		normals[2][1],
		normals[2][2],
		1
	};
	applyTransformation(Xnorm[matlevel - 1], homogenousN0, &resultX, &resultY, &resultZ, &resultW);
	transformedNormals[0][0] = resultX / resultW;
	transformedNormals[0][1] = resultY / resultW;
	transformedNormals[0][2] = resultZ / resultW;
	//normalise
	float length = calculateNorm(transformedNormals[0]);
	transformedNormals[0][0] = transformedNormals[0][0] / length;
	transformedNormals[0][1] = transformedNormals[0][1] / length;
	transformedNormals[0][2] = transformedNormals[0][2] / length;

	applyTransformation(Xnorm[matlevel - 1], homogenousN1, &resultX, &resultY, &resultZ, &resultW);
	transformedNormals[1][0] = resultX / resultW;
	transformedNormals[1][1] = resultY / resultW;
	transformedNormals[1][2] = resultZ / resultW;
	//normalise
	length = calculateNorm(transformedNormals[0]);
	transformedNormals[1][0] = transformedNormals[1][0] / length;
	transformedNormals[1][1] = transformedNormals[1][1] / length;
	transformedNormals[1][2] = transformedNormals[1][2] / length;

	applyTransformation(Xnorm[matlevel - 1], homogenousN2, &resultX, &resultY, &resultZ, &resultW);
	transformedNormals[2][0] = resultX / resultW;
	transformedNormals[2][1] = resultY / resultW;
	transformedNormals[2][2] = resultZ / resultW;
	//normalise
	length = calculateNorm(transformedNormals[0]);
	transformedNormals[2][0] = transformedNormals[2][0] / length;
	transformedNormals[2][1] = transformedNormals[2][1] / length;
	transformedNormals[2][2] = transformedNormals[2][2] / length;

	/**************************** NEW ADDITION HW6 *************************************************/
	for (int i = 0; i < 3; i++)
		transformedVertices[i][0] = transformedVertices[i][0] - offsetX[rendNo];
	for (int i = 0; i < 3; i++)
		transformedVertices[i][1] = transformedVertices[i][1] - offsetY[rendNo];
	
	rasterize(transformedVertices, transformedNormals, texels);

	return GZ_SUCCESS;
}

