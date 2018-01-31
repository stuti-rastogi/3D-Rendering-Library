#include	"gz.h"
#ifndef GZRENDER_
#define GZRENDER_


/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

class GzRender {			/* define a renderer */


public:
	unsigned short	xres;
	unsigned short	yres;
	GzPixel		*pixelbuffer;		/* frame buffer array */
	char* framebuffer;

	GzCamera		m_camera;
	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
	GzColor		flatcolor;          /* color state for flat shaded triangles */
	int			interp_mode;
	int			numlights;
	GzLight		lights[MAX_LIGHTS];
	GzLight		ambientlight;
	GzColor		Ka, Kd, Ks;
	float		    spec;		/* specular power */
	GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */

/**************************** NEW ADDITION HW6 *************************************************/
	float		offsetX[AAKERNEL_SIZE];
	float		offsetY[AAKERNEL_SIZE];
	//float		weight[AAKERNEL_SIZE];
	int			rendIndex;

								// Constructors
	GzRender(int xRes, int yRes);
	~GzRender();

	// Function declaration

	// HW1: Display methods
	int GzDefault();
	int GzBeginRender();
	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth	*z);

	int GzFlushDisplay2File(FILE* outfile);
	int GzFlushDisplay2FrameBuffer();

	// HW2: Render methods
	int GzPutAttribute(int numAttributes, GzToken *nameList, GzPointer *valueList);
	int GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList, int rendNo);

	// HW3
	int GzPutCamera(GzCamera camera);
	int GzPushMatrix(GzMatrix	matrix);
	int GzPopMatrix();

	// Extra methods: NOT part of API - just for general assistance */
	inline int ARRAY(int x, int y) { return (x + y*xres); }	/* simplify fbuf indexing */
	inline short	ctoi(float color) { return(short)((int)(color * ((1 << 12) - 1))); }		/* convert float color to GzIntensity short */


																								// Object Translation
	int GzRotXMat(float degree, GzMatrix mat);
	int GzRotYMat(float degree, GzMatrix mat);
	int GzRotZMat(float degree, GzMatrix mat);
	int GzTrxMat(GzCoord translate, GzMatrix mat);
	int GzScaleMat(GzCoord scale, GzMatrix mat);

private:
	bool isHorizontal(float y0, float y1);
	int assignFlavorsHorizontal(GzCoord vertices[], int hor_v1, int hor_v2, int other,
		int *e01_flavor, int *e12_flavor, int *e20_flavor);
	int assignFlavorsGeneral(GzCoord vertices[], float ymin, float ymax,
		int *e01_flavor, int *e12_flavor, int *e20_flavor);
	void computeBoundingBox(float x0, float x1, float x2, float y0, float y1, float y2,
		float *xmin, float *xmax, float *ymin, float *ymax);
	void computeEquation(float tailX, float tailY, float headX, float headY,
		float *a, float *b, float *c);
	void getCWEdgesHori(GzCoord vertices[], int hori_v1, int hori_v2, int other, int hori_flavor,
		int *e1_head, int *e1_tail, int *e2_head, int *e2_tail, int *e3_head, int *e3_tail);
	void getCWEdgesGeneral(GzCoord vertices[], int soleEdge_flavor,
		int minV, int maxV, int other, int *e1_head, int *e1_tail, int *e1_flavor,
		int *e2_head, int *e2_tail, int *e2_flavor, int *e3_head, int *e3_tail, int *e3_flavor);
	int evaluateEquation(float a, float b, float c, float x, float y);

	//hw3 helper methods
	float calculateNorm(GzCoord vector);
	void calculateCrossProduct(GzCoord a, GzCoord b, float *resultX, float *resultY, float *resultZ);
	float calculateDotProduct(GzCoord a, GzCoord b);
	void multiplyMatrices(GzMatrix mat1, GzMatrix mat2, GzMatrix * result);
	void applyTransformation(GzMatrix stack, float vertex[4],
		float* resultx, float *resulty, float *resultz, float *resultw);
	void rasterize(GzCoord transformedVertices[], GzCoord transformedNormals[], GzTextureIndex texels[]);

	//hw4 additions
	int GzPushMatrixOnNormStack(GzMatrix matrix);
	void scalarVectorMultiply(float a, GzCoord v, float *resultX, float *resultY, float *resultZ);
	void calculateDifferenceVectors(GzCoord a, GzCoord b, float *resultX, float *resultY, float *resultZ);
	void computeShadingEquation(GzCoord transformedNormal, float *r, float *g, float *b,
		GzColor Ka_tex, GzColor Kd_tex, GzColor Ks_tex);
	int checkNormalSide(GzCoord normal, GzCoord light, GzCoord E, float *nx, float *ny, float *nz);
	float interpolate(GzCoord v1, GzCoord v2, GzCoord v3, int x, int y);

};
#endif