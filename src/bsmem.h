/* Header for oifits routines.
 * This package uses the Oifits Exchange routines by John Young to view,
 * select and extract oi data.
 *
 * Hrobjartur Thorsteinsson 9/12/03
 * Fabien Baron 2004-2009
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <termios.h>
#include <complex.h>
#include <fftw3.h>
#include <nfft3.h>
#include "exchange.h"
#include "fitsio.h"
//#include "random.h"



#define maxins 50
#define billion (1.0e9)

/**********************/
/* User's application */
/**********************/

/*
  Structure for carrying information
  on data format and analysis
*/

typedef struct 
{
  // OIFITS SELECTION
  char datafile[FLEN_FILENAME];
  int target_id;
  char targetname[FLEN_VALUE];
  long numvis2;
  long numt3;
  double wavel;
  char insname[maxins][FLEN_VALUE];
  long numins;
  float minband;
  float maxband;
  // RECONST PARAMETERS  
  float modelwidth;
  int modeltype;
  float modelflux;
  float blur;
  float alpha;
  int t3io;
  int v2io;
  int phsio;
  int box;
  int xydim;
  float xyint;
  int maxiter;
  int entropy;
  int regularization;
  int verbose;
  int noui;
  char fitsimage[100];
  char priorfile[100];
  float displaypower;
  int ntelescopes;
  float fluxerr;
  int t3errtype;
  int forced_extrapolation;
}RECONST_PARAMETERS;


typedef struct
{  /* Space */
  int iNv2;
  int iNt3;
  int iNUV;
  int iNX;
  int iNY;
  int DFT;
  float xyint;
  float fblur;
  /* Work environments */
  float complex *cvis;
  //  float complex *acT0;
  float complex *dcvis;
  float complex *data_phasor;
  // float complex *DFT_table;
  float imflux;
  float difflux;
  float *uvpoint;
  int idimage;
  int idv2;
  int idt3phi;
  int idt3amp;
  // int *aiTRI;
  nfft_plan p;
}USER;


#define UINFO_FAIL          -12
#define D_OPEN_ERR			-20
#define D_READ_ERR			-21
#define D_WRIT_ERR			-22
#define M_ALLOC_ERR			-30
int s_equals(char *A, char *B);

int get_default_model(float* model, float* currentimage, int n, float xs, float xe, float ys, float ye, float xyint, int id);
int get_box(float* xr,float* yr,float* x,float* y,int id);
int redisplay(float* image, USER *user, float displaypower, int contour);
int open_redisplay(USER *user, char* dev);
int close_redisplay(USER *user );
int dispuv(char* dev);
void center_img( float* img,int N );
float abs2(float complex cp);

/* Function declarations */
int get_oi_fits_data(RECONST_PARAMETERS *reconst_parameters, int* status);
int read_fits_image(char* fname, float* img, int* n, float* xyint, char* source, char* datafile, int* status);

int import_single_epoch_oifits( char *filename, bool use_v2, bool use_t3amp, bool use_t3phi, bool use_visamp, bool use_visphi,
                               double v2a, double v2s, double t3ampa, double t3amps, double t3phia, double t3phis,
                               double visampa, double visamps, double visphia, double visphis, double fluxs, double cwhm,
				double uvtol, int nwavr, double *wavmin, double *wavmax, double *timemin, double *timemax);

/* Error exit handling */
#define CALL(x) {if((CALLvalue = (x))<0) return CALLvalue;}
#define SUCCESS 0
/* Constants */
#define M_PI 3.14159265358979323
#define MAS (3.14159265358979323/180.0)/3600000.0
#define OVERSAMPLING 3.0 /* compared to Shannon. 6 pix/fringe -> 3.0 */
#define infinity 1.0e9
#ifndef MSVS 
        #define FTEST ftest_
        #define UVCHK uvchk_
        #define OPSET opset_
        #define OP1 op1_
        #define TROP1 trop1_
	#define UAREA uarea_
	#define MEINIT meinit_
	#define MEMSET memset_
	#define MEMGET memget_
	#define USETM usetm_
	#define UGETM ugetm_
	#define MEREST merest_
	#define MESAVE mesave_
	#define MEMTRQ memtrq_
	#define MEM4 mem4_
	#define MOVIE4 movie4_
        #define MASK4 mask4_
        #define MEMEX memex_
        #define VMEMEX vmemex_
	#define VOPUS vopus_
	#define VTROP vtrop_
	#define USAVE usave_
	#define UREST urest_
	#ifdef __stdcall
		#undef __stdcall
	#endif
	#define __stdcall 
#endif
