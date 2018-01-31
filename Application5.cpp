// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "ppot.asc"
#define OUTFILE "output.ppm"

extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */
extern int GzFreeTexture();

void shade(GzCoord norm, GzCoord color);

/**************************** NEW ADDITION HW6 *************************************************/
float AAFilter[AAKERNEL_SIZE][3] = { -0.52, 0.38, 0.128,                  0.41, 0.56, 0.119,                     0.27, 0.08, 0.294,
-0.17, -0.29, 0.249,                    0.58, -0.55, 0.104,                   -0.31, -0.71, 0.106 };

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = 512;		// frame buffer and display width
	m_nHeight = 512;    // frame buffer and display height

	m_pRender = new GzRender(m_nWidth, m_nHeight);
	m_pRender->GzDefault();

	m_pFrameBuffer = m_pRender->framebuffer;

		//creating 6 different renderers
		/*renderArray[rendNo] = new GzRender(m_nWidth, m_nHeight);
		renderArray[rendNo]->GzDefault();

		frameBufferArray[rendNo] = renderArray[rendNo]->framebuffer;*/

		/* Translation matrix */
		GzMatrix	scale =
		{
			3.25,	0.0,	0.0,	0.0,
			0.0,	3.25,	0.0,	-3.25,
			0.0,	0.0,	3.25,	3.5,
			0.0,	0.0,	0.0,	1.0
		};

		GzMatrix	rotateX =
		{
			1.0,	0.0,	0.0,	0.0,
			0.0,	.7071,	.7071,	0.0,
			0.0,	-.7071,	.7071,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

		GzMatrix	rotateY =
		{
			.866,	0.0,	-0.5,	0.0,
			0.0,	1.0,	0.0,	0.0,
			0.5,	0.0,	.866,	0.0,
			0.0,	0.0,	0.0,	1.0
		};

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
		camera.position[X] = -3;
		camera.position[Y] = -25;
		camera.position[Z] = -4;

		camera.lookat[X] = 7.8;
		camera.lookat[Y] = 0.7;
		camera.lookat[Z] = 6.5;

		camera.worldup[X] = -0.2;
		camera.worldup[Y] = 1.0;
		camera.worldup[Z] = 0.0;

		camera.FOV = 63.7;              /* degrees *              /* degrees */

		//status |= renderArray[rendNo]->GzPutCamera(camera);
		status |= m_pRender->GzPutCamera(camera);
#endif 

		/* Start Renderer */
		//status |= renderArray[rendNo]->GzBeginRender();
		status |= m_pRender->GzBeginRender();

		/* Light */
		GzLight	light1 = { { -0.7071, 0.7071, 0 },{ 0.5, 0.5, 0.9 } };
		GzLight	light2 = { { 0, -0.7071, -0.7071 },{ 0.9, 0.2, 0.3 } };
		GzLight	light3 = { { 0.7071, 0.0, -0.7071 },{ 0.2, 0.7, 0.3 } };
		GzLight	ambientlight = { { 0, 0, 0 },{ 0.3, 0.3, 0.3 } };

		/* Material property */
		GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
		GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
		GzColor diffuseCoefficient = { 0.7, 0.7, 0.7 };

		/*
		renderer is ready for frame --- define lights and shader at start of frame
		*/

		/*
		* Tokens associated with light parameters
		*/
		nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[0] = (GzPointer)&light1;
		nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[1] = (GzPointer)&light2;
		nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
		valueListLights[2] = (GzPointer)&light3;
		//status |= renderArray[rendNo]->GzPutAttribute(3, nameListLights, valueListLights);
		status |= m_pRender->GzPutAttribute(3, nameListLights, valueListLights);

		nameListLights[0] = GZ_AMBIENT_LIGHT;
		valueListLights[0] = (GzPointer)&ambientlight;
		//status |= renderArray[rendNo]->GzPutAttribute(1, nameListLights, valueListLights);
		status |= m_pRender->GzPutAttribute(1, nameListLights, valueListLights);

		/*
		* Tokens associated with shading
		*/
		nameListShader[0] = GZ_DIFFUSE_COEFFICIENT;
		valueListShader[0] = (GzPointer)diffuseCoefficient;

		/*
		* Select either GZ_COLOR or GZ_NORMALS as interpolation mode
		*/

		for (int rendNo = 0; rendNo < AAKERNEL_SIZE; rendNo++)
		{
			nameListShader[1] = GZ_INTERPOLATE;
			//interpStyle = GZ_COLOR;         /* Gouraud shading */
			interpStyle = GZ_NORMALS;         /* Phong shading */
			valueListShader[1] = (GzPointer)&interpStyle;

			nameListShader[2] = GZ_AMBIENT_COEFFICIENT;
			valueListShader[2] = (GzPointer)ambientCoefficient;
			nameListShader[3] = GZ_SPECULAR_COEFFICIENT;
			valueListShader[3] = (GzPointer)specularCoefficient;
			nameListShader[4] = GZ_DISTRIBUTION_COEFFICIENT;
			specpower = 32;
			valueListShader[4] = (GzPointer)&specpower;

			nameListShader[5] = GZ_TEXTURE_MAP;
#if 0   /* set up null texture function or valid pointer */
			valueListShader[5] = (GzPointer)0;
#else
			valueListShader[5] = (GzPointer)(tex_fun);	/* or use ptex_fun */
#endif

		/**************************** NEW ADDITION HW6 *************************************************/
			nameListShader[6] = GZ_REND_NO;
			valueListShader[6] = (GzPointer)&rendNo;
			nameListShader[7] = GZ_AASHIFTX;
			valueListShader[7] = (GzPointer)&AAFilter[rendNo][0];
			nameListShader[8] = GZ_AASHIFTY;
			valueListShader[8] = (GzPointer)&AAFilter[rendNo][1];
			//status |= renderArray[rendNo]->GzPutAttribute(9, nameListShader, valueListShader);
			status |= m_pRender->GzPutAttribute(9, nameListShader, valueListShader);
		}


		//status |= renderArray[rendNo]->GzPushMatrix(scale);
		//status |= renderArray[rendNo]->GzPushMatrix(rotateY);
		//status |= renderArray[rendNo]->GzPushMatrix(rotateX);
		status |= m_pRender->GzPushMatrix(scale);
		status |= m_pRender->GzPushMatrix(rotateY);
		status |= m_pRender->GzPushMatrix(rotateX);

		if (status) exit(GZ_FAILURE);

		if (status)
			return(GZ_FAILURE);
		else
			return(GZ_SUCCESS);
}

int Application5::Render()
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */
	GzCoord		normalList[3];	/* vertex normals */
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */
	char		dummy[256];
	int			status;

	int resolution = m_nHeight * m_nWidth;
	//new pixel buffer to accumulate filtered values in application
	GzPixel *pbFilter = new GzPixel[resolution];

	//set all new pixel buffer values to 0
	for (int i = 0; i < resolution; i++)
	{
		pbFilter[i].red = 0;
		pbFilter[i].green = 0;
		pbFilter[i].blue = 0;
		pbFilter[i].alpha = 0;
		pbFilter[i].z = 0;
	}

	//declare variables for GzGet
	GzIntensity r, g, b, a;
	GzDepth z;

	/* Initialize Display */
	status |= m_pRender->GzDefault();  /* init for new frame */

	for (int rendNo = 0; rendNo < AAKERNEL_SIZE; rendNo++)
	{
		//status |= renderArray[rendNo]->GzDefault();  /* init for new frame */
		status |= m_pRender->GzDefault();
		/*
		* Tokens associated with triangle vertex values
		*/
		nameListTriangle[0] = GZ_POSITION;
		nameListTriangle[1] = GZ_NORMAL;
		nameListTriangle[2] = GZ_TEXTURE_INDEX;

		// I/O File open
		FILE *infile;
		if ((infile = fopen(INFILE, "r")) == NULL)
		{
			AfxMessageBox("The input file was not opened\n");
			return GZ_FAILURE;
		}

		/*
		* Walk through the list of triangles, set color
		* and render each triangle
		*/
		while (fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[0][0]), &(vertexList[0][1]),
				&(vertexList[0][2]),
				&(normalList[0][0]), &(normalList[0][1]),
				&(normalList[0][2]),
				&(uvList[0][0]), &(uvList[0][1]));
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[1][0]), &(vertexList[1][1]),
				&(vertexList[1][2]),
				&(normalList[1][0]), &(normalList[1][1]),
				&(normalList[1][2]),
				&(uvList[1][0]), &(uvList[1][1]));
			fscanf(infile, "%f %f %f %f %f %f %f %f",
				&(vertexList[2][0]), &(vertexList[2][1]),
				&(vertexList[2][2]),
				&(normalList[2][0]), &(normalList[2][1]),
				&(normalList[2][2]),
				&(uvList[2][0]), &(uvList[2][1]));

			/*
			 * Set the value pointers to the first vertex of the
			 * triangle, then feed it to the renderer
			 * NOTE: this sequence matches the nameList token sequence
			 */
			valueListTriangle[0] = (GzPointer)vertexList;
			valueListTriangle[1] = (GzPointer)normalList;
			valueListTriangle[2] = (GzPointer)uvList;
			m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle, rendNo);
		}

		//accumulate the new renderer pixel buffer values in application pixel buffer
		float weight = AAFilter[rendNo][2];

		for (int y = 0; y < m_nHeight; y++)
		{
			for (int x = 0; x < m_nWidth; x++)
			{
				m_pRender->GzGet(x, y, &r, &g, &b, &a, &z);
				pbFilter[y*m_nWidth + x].red = pbFilter[y*m_nWidth + x].red + (weight * r);
				pbFilter[y*m_nWidth + x].green = pbFilter[y*m_nWidth + x].green + (weight * g);
				pbFilter[y*m_nWidth + x].blue = pbFilter[y*m_nWidth + x].blue + (weight * b);
				pbFilter[y*m_nWidth + x].alpha = pbFilter[y*m_nWidth + x].alpha + (weight * a);
				pbFilter[y*m_nWidth + x].z = pbFilter[y*m_nWidth + x].z + (weight * z);
			}
		}

		if (fclose(infile))
			AfxMessageBox(_T("The input file was not closed\n"));
	}	//end of for loop for subsamples

	//now put accumulated pixel buffer back into renderer pixel buffer for flushing
	for (int y = 0; y < m_nHeight; y++)
	{
		for (int x = 0; x < m_nWidth; x++)
		{
			GzPixel current = pbFilter[y*m_nWidth + x];
			m_pRender->GzPut(x, y, current.red, current.green, current.blue, current.alpha, current.z);
		}
	}

	FILE *outfile;
	if ((outfile = fopen(OUTFILE, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
		return GZ_FAILURE;
	}

	m_pRender->GzFlushDisplay2File(outfile); 	/* write out or update display to file*/
	m_pRender->GzFlushDisplay2FrameBuffer();	// write out or update display to frame buffer

	/* 
	 * Close file
	 */ 

	if( fclose( outfile ) )
      AfxMessageBox(_T( "The output file was not closed\n" ));
 
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	free(m_pRender);
	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



