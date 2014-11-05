/* Header for oifits routines.
 * This package uses the Oifits Exchange routines by John Young to view,
 * select and extract oi data.
 *
 * Hrobjartur Thorsteinsson 9/12/03
 * Fabien Baron 2004-2009
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <complex.h>
#include <fftw3.h>
#include <nfft3.h>
#include "exchange.h"
#include "fitsio.h"
#include "random.h"

#define maxins 50
#define billion (1.0e9)

typedef struct _uvpnt
{
  short sign;
  long uvpnt;
}oi_uvpnt;

typedef struct _bsref
{
  /* Structure for bisp to uv coord table referencing.
   * Negative uv table number means the particular baseline is conjugated.
   */
  oi_uvpnt ab;
  oi_uvpnt bc;
  oi_uvpnt ca;
}oi_bsref;

typedef struct _uv
{
  float u;
  float v;
  float bandwidth;
  float wavelength;
}oi_uv;

typedef struct _data
{
  float *cvis;
  float *pow;
  float *powerr;
  float *bisamp;
  float *bisamperr;
  float *bisphs;
  float *bisphserr;
  oi_uv *uv;
  oi_bsref *bsref;
  int npow;
  int nbis;
  int nuv;
}oi_data;

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
  int bisio;
  int powio;
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
  float v2a;
  float v2b;
  float t3ampa;
  float t3ampb;
  float t3phia;
  float t3phib;
  int novis2;
  int not3amp;
  int not3phi;
  float displaypower;
  int ntelescopes;
  float fluxerr;
  int biserrtype;
  int forced_extrapolation;
}RECONST_PARAMETERS;

typedef struct {
  float u;
  float v;
}AmpUV;

typedef struct {
  float u1;
  float v1;
  float u2;
  float v2;
}BsUV;

typedef struct
{  /* Space */
  int iNpow;
  int iNbis;
  int iNUV;
  int iNX;
  int iNY;
  int DFT;
  float xyint;
  float fblur;
  /* Work environments */
  float complex *current_visi;
  float complex *acT0;
  float complex *current_diffvisi;
  float complex *data_phasor;
  // float complex *DFT_table;
  float imflux;
  float difflux;
  int *filter_powerspectrum;
  int *filter_bispectrum;
  oi_bsref *bsref;
  /* Work environment for transform/gridding routines */
  float *UV_point;
  float *filter_UV;
  int idimage;
  int idpow;
  int idt3phs;
  int idt3amp;
  // int *aiTRI;
  nfft_plan p;
}USER;


#define UINFO_FAIL          -12
#define D_OPEN_ERR			-20
#define D_READ_ERR			-21
#define D_WRIT_ERR			-22
#define M_ALLOC_ERR			-30

int get_default_model(float* model, float* currentimage, int n, float xs, float xe, float ys, float ye, float xyint, int id);
int get_box(float* xr,float* yr,float* x,float* y,int id);
int redisplay(float* image, USER *user, oi_data *data, float displaypower, int contour);
int open_redisplay(USER *user, char* dev);
//int close_redisplay(int id);
int dispuv(oi_uv* A, oi_bsref* B, int nuv, int npow, int nbis, char* dev);
void center_img( float* img,int N );
float abs2(float complex cp);

/* Function declarations */
int get_oi_fits_selection(RECONST_PARAMETERS *reconst_parameters, int* status);
int get_oi_fits_data(RECONST_PARAMETERS *reconst_parameters, oi_data *data, int* status);
int compare_uv(oi_uv uv, oi_uv withuv, float thresh);
void free_oi_target(oi_target *targets);
void free_oi_wavelength(oi_wavelength *wave);
void free_oi_vis2(oi_vis2 *vis2);
void free_oi_t3(oi_t3 *t3);
void free_oi_data(oi_data *data);
int count_redundant_bsuv(oi_bsref *bsref, int nbs);
float bsuv_coverage_quality(oi_bsref *bsref, int nbs, oi_uv *uv, int nuv);
int read_fits_image(char* fname, float* img, int* n, float* xyint, char* source, char* datafile, int* status);
int close_redisplay(USER *user );

/* Error exit handling */
#define CALL(x) {if((CALLvalue = (x))<0) return CALLvalue;}
#define SUCCESS 0
/* Constants */
#define pi 3.14159265358979323
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
