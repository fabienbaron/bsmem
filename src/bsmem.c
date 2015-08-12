/* bsmem.c
 *
 * The Contents of this file are made available subject to the terms
 * of the GNU Lesser General Public License:
 *
 * Copyright ( C ) 2002-2015 Fabien Baron, David Buscher, Hrobjartur Thorsteinsson
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or ( at your option ) any later version.
 *
 * This file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file. If not, see
 * http://www.gnu.org/licenses/
 */

/* * Project:   BSMEM
 *
 *
 * Purpose:   To invert optical aperture synthesis data consisting of visiblity
 *       amplitudes and triple products ( bispectrum measurements )
 *
 * History:   This is the original implementation of BSMEM created by David F. Buscher
 *       in Fortran in 1990, but rewritten ( with permission ) by
 *       Hrobjartur Thorsteinsson using the new memsys4 library in F.
 *       Fabien Baron took over development in mid 2004.
 *       C to Fortran translation and toy program provided by Steve Gull.
 */

/*
 * STATUS: 01/02/03 Program working.
 *     09/12/03 Added oi-fits routines for reading oi-fits formated data.
 *     15/12/03 Pixel weightings == model. MEMSYS defined.
 *     19/04/06 Everything compiles without warnings under -Wall.
 *     01/09/07 Numerous bugs fixed. Added the ERRMAP command to create an error map.
 *     10/09/07 Move to <complex.h> C library ( much faster ).
 *     01/02/08 Many minor I/O bugs fixed.
 *     01/05/08 Further minor fixes for Beauty Contest, help corrected.
 *     12/07/08 Fixes for GCC4 + ISO 99 compliance.
 *     01/02/09 Better default bispectrum noise model. Gridding now done with NFFT. More prior options.
 *     01/06/09 Extrapolation of powerspectra to create pseudo-data on triple amolitudes. v.1.3
 *     10/01/10 Forced normalization by introducing a zero flux powerspectrum. v1.4
 *     20/03/10 Added triple amplitude visu, minor fixes, waveband selection from commandline. v1.5
 *     13/09/11 Fixed memory leaks using cppcheck, fixed full elliptic approximation of bispectrum errors v1.6
 *     10/11/14 Updated convexification based on WISARD equations
 *     05/08/15 Minor tweaks for PFI
 *     11/08/15 Major revision: OIFITS2 loader from SQUEEZE, complex VIS support v2.0
 */

#include "bsmem.h"
#include "cpgplot.h"

int oi_hush_errors = 0; // flag for read_fits.c

extern void __stdcall UAREA( int* , int* , int* , int* , int* , int* );
extern void __stdcall MEINIT( );
extern void __stdcall MEMSET( int* , int* , int* , float* , float* , float* , float* , float* );
extern void __stdcall MEMGET( int* , int* , int* , float* , float* , float* , float* , float* );
extern void __stdcall USETM( float* , int* , float* , int* );
extern void __stdcall UGETM( float* , int* , float* , int* );
extern void __stdcall MEMEX( float* , int* , int* );
extern void __stdcall MEMTRQ( float* , int* , float* );
extern void __stdcall MEREST( );
extern void __stdcall MESAVE( );
extern void __stdcall MEM4( float* , int* , float* , float* , float* , float* , float* , float* , float* , float* , float* , float* ,
    float* , float* , int* , int* );
extern void __stdcall MOVIE4( float* , int* , int* , int* );
extern void __stdcall MASK4( float* , float* , float* , float* , int* , int* );
void __stdcall VTROP( float *dh , float *dd );
void __stdcall VOPUS( float *dh , float *dd );
void __stdcall VMEMEX( float *hidden , float *data );

// Subroutines
int iterate( int niters , RECONST_PARAMETERS* reconst_parameters );
float square( float n );
int read_model( char* filename , float* mod , int N , float normalisation );
int read_data( FILE* f1 );
int get_input( RECONST_PARAMETERS* reconst_parameters , int argc , char** argv );
int compare( char* , char* );
int bsmem( int argc , char** argv );
int set_memsys_dataspace( float *st , int *kb ,  RECONST_PARAMETERS *reconst_parameters );
int reset_memsys_dataspace( float *st , int *kb ,  RECONST_PARAMETERS *reconst_parameters );
int write_fits_image( float* img , RECONST_PARAMETERS* reconst_parameters , int* status );
void uvchuck(float *st , int *kb , RECONST_PARAMETERS *reconst_parameters );
void display_oifits_info( );
void set_model( RECONST_PARAMETERS* reconst_parameters , float *startmod );
void display_command_help( );
void op_trop_check( int npix , int ndata );

// DFT precomputed coefficient table
float complex* DFT_table = NULL;

// Global variables
static float *st;
int kb[ 40 ];
// unsigned Rand[ 113 ];
USER user;
int nscales = 4;

float plow, phigh,alpha ;
int total_iter;

/* Global OIFITS variables */
long nuv;
char *oifits_file;
long nvis, nv2, nt3, nt3phi, nt3amp, nt3amp_orphans, nt3phi_orphans, nvisamp, nvisphi, nvisamp_orphans, nvisphi_orphans;
long *visin, *v2in, *t3in1, *t3in2, *t3in3;
double *u, *v;
int ntimer;
bool use_diffvis = FALSE; // default for VIS tables = complex vis, not differential vis
long *dvisnwav;
long **dvisindx;
double *__restrict uv_time;
double *__restrict uv_lambda;
double *__restrict uv_dlambda;
int *uvwav2chan = NULL;
int *uvtime2chan = NULL;
double *__restrict data;
double *__restrict data_err;

int main( int argc , char** argv )
{
  int callvalue;
  callvalue = bsmem(argc, argv);
  if (callvalue != 0)
    printf("\n\nERROR CODE = %d\n", callvalue);
  return SUCCESS;
}
float sinc( float x)
{
  if( x < 1e-7)
    return 1.;
  else
    return sin(x)/x;
}

int bsmem( int argc , char** argv )
{
  int memrun, samples, istat, kstat, ntrnsx, ncorr, ntrans, level = 10, method[ 4 ], nrand = 1, iseed = 1234;
  float def = -1.0, aim = 1.0, rate = 1.8, acc = -1.0, utol = 0.1;
  int status;
  int disp_id = 0;
  float entropy, test, chisq, omega;
  float scale;
  float pdev, glow, ghigh, gdev;

  /* Loop variables */
  int i, j, ii;//, jj, uu;
  int stop;
  int ndata, ncells, nstore;
  float *startmod, *maps, *errmap, *meanmap, *temp;
  float sum;

  /* Input structure */
  RECONST_PARAMETERS reconst_parameters;

  /* Parameters for MEMSYS user interface */
  char line[ 150 ];
  char sc1[ 50 ], sc2[ 50 ];
  int si2;
  char tempstr;
  int contour = 1;
  /* OI FITS data structure */
  status = 0;

  /* Display variables */
  char dispdev[ 30 ];

  /**************** Get user input parameters for reconstruction *******/
  if (get_input(&reconst_parameters, argc, argv) != 0)
    goto EXIT;

  /**************** Graphic display initialization ********************/
  strcpy(dispdev, "/xwindow");
  reconst_parameters.displaypower = 1.0;

  /**************** Read OIFITS formated file *************************/
  printf("Datafile:\t\t%s\n", reconst_parameters.datafile);
  get_oi_fits_data(&reconst_parameters, &status);

  /**************** Map dimensions and parameters ************************************/
  user.iNX = reconst_parameters.xydim;
  user.iNY = reconst_parameters.xydim;
  user.xyint = reconst_parameters.xyint;

  /**************** Allocate MEMSYS storage areas *************************************/
  ndata = 2 * nvis + nv2 + 2 * nt3;
  ncells = user.iNX * user.iNY;
  nstore = 6 * ncells + 9 * ndata;
  UAREA(&ncells, &ncells, &ndata, &nstore, &level, kb);
  st = malloc(nstore * sizeof(float));

  /*************** Assign memory to work spaces *********************************/
  /* Opus/Tropus */
  user.current_cvis = malloc( nuv * sizeof( float complex ) );
  user.current_derivcvis = malloc( nuv * sizeof( float complex ) );
  user.data_phasor = malloc( nt3 * sizeof( float complex ) );
  
  for (i = 0; i < nuv; i++)
  {
    user.current_cvis[ i ] = 0.0;
    user.current_derivcvis[ i ] = 0.0;
  }

  
  /*************** Initialise UV_point and t3in for transforms *********/
  
  printf("Number of unique uv points in NFFT: %ld\n", nuv);

  /* Convert exp( -I*closure phase ) to Real/Imaginary parts */
  for (ii = 0; ii < nt3; ii++)
    user.data_phasor[ ii ] = cexp(-I * data[ nv2 + nt3amp + ii ] * M_PI / 180.0);

  /* Set up MEMSYS data area*/
  set_memsys_dataspace(st, kb,  &reconst_parameters);
  uvchuck(st, kb, &reconst_parameters);

  int NN[ 2 ], nn[ 2 ];
  NN[ 0 ] = user.iNX;
  nn[ 0 ] = 2 * NN[ 0 ];
  NN[ 1 ] = user.iNX;
  nn[ 1 ] = 2 * NN[ 1 ];
  nfft_init_guru(&user.p, 2, NN, nuv, nn, 6, PRE_FULL_PSI | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  for (i = 0; i < nuv; i++)
  {
    user.p.x[ 2 * i ] = v[ i ] * user.xyint * MAS;
    user.p.x[ 2 * i + 1 ] = u[ i ] * user.xyint * MAS;
  }

  nfft_precompute_one_psi(&user.p);
  
  /*************** Compulsory initialisation ***********************************/
  MEINIT();

  if (reconst_parameters.verbose == 1)
  {
    printf("Entropy functional:\t");
    if (reconst_parameters.entropy == 1)
      printf("Gull-Skilling entropy\n");
     if (reconst_parameters.entropy == 2)
      printf("Positive/negative Gull-Skilling entropy\n");
    if (reconst_parameters.entropy == 4)
      printf("Quadratic/L2\n");
    printf("Hyperparameter scheme:\t");
    if (reconst_parameters.regularization == 1)
      printf("Classic Bayesian\n");
    if (reconst_parameters.regularization == 2)
      printf("Classic Bayesian with noise scaling\n");
    if (reconst_parameters.regularization == 3)
      printf("Fixed alpha=%f\n", reconst_parameters.alpha);
    if (reconst_parameters.regularization == 4)
      printf("Chi2 = N method\n");
    printf("Maximum n# iterations:\t");
    if (reconst_parameters.maxiter <= 0)
      printf("Infinite\n");
    else
      printf("%d\n", reconst_parameters.maxiter);
  }

  method[ 0 ] = reconst_parameters.regularization;
  method[ 1 ] = reconst_parameters.entropy;
  method[ 2 ] = 1; /* Gaussian noise */
  method[ 3 ] = 2; /* for non-linearity */

  def = -1.0;
  if (reconst_parameters.regularization == 3)
  {
    aim = reconst_parameters.alpha;
  }
  else
    aim = 1.0;

  rate = 1.;

  MEMSET(method, &nrand, &iseed, &aim, &rate, &def, &acc, &utol);

  // Set up model
  startmod = malloc(user.iNX * user.iNY * sizeof(float));
  set_model(&reconst_parameters, startmod);

  RESTART:

  /**************** Set pixel-weighting ***************************************/
  for (i = 0; i < user.iNX * user.iNY; i++)
    st[ kb[ 2 ] + i ] = startmod[ i ];

  for (i = 0; i < user.iNX * user.iNY; i++)
    st[ kb[ 0 ] + i ] = startmod[ i ]; /*kb[ 0 ] points at FORTRAN area ST( 1 )*/

  /**************** Set V0, VIS etc... ***************************************************/
  temp = malloc(ndata * sizeof(float));
  VMEMEX(startmod, temp);
  free(temp);

   /**************** MaxEnt iterations *******************************************/
  if (reconst_parameters.verbose == 1)
    printf("\nStarting Maximum Entropy Reconstruction.\n");
  memrun = 1;
  total_iter = 1;
  MEM4(st, &memrun, &entropy, &test, &chisq, &scale, &plow, &phigh, &pdev, &glow, &ghigh, &gdev, &omega, &alpha, &istat, &ntrans);
  if (reconst_parameters.verbose == 1)
  {
    printf("Iteration %d Ntrans === %d istat === %d%d%d%d%d%d%d \n", total_iter, ntrans, istat % 2, istat / 2 % 2, istat / 4 % 2, istat / 8
        % 2, istat / 16 % 2, istat / 32 % 2, istat / 64 % 2);
    printf("Entropy === %f  Chisq === %f Flux === %f Alpha === %f Omega === %f \n", entropy, chisq / (float) ndata, user.imflux, alpha, omega);
    printf("Logprob === %f Good Measurements === %f Scale === %f \n", (plow + phigh) / 2., (glow + ghigh) / 2., scale);
  }
  memrun = 2;

  int ninterface_commands = 17;
  char interface_commands[ ][ 10 ] =
  { "EXIT", "CENTER", "SNR", "RESTART", "UV", "HELP", "STARTMOD", "ERRMAP", "MEANMAP", "WRITEFITS", "SCALE", "READMOD", "VERBOSE",
      "REDISP", "DO", "CONTOUR", "MISSING" };

  int command_number;
  int exit = 0;
  int n;

  // Command selection
  do
  {
    command_number = ninterface_commands + 1;
    if (reconst_parameters.noui == 0)
    {
      printf("* ");
      if( fgets(line, 150, stdin) == NULL)
	printf("Error getting command\n");
      tempstr = line[ 0 ];
      n = sscanf(line, "%s %s", sc1, sc2);
      if (strlen(line) == (strlen(sc1) + 1))
        strcpy(sc2, "");
      for (i = 0; i < ninterface_commands; i++)
        if (strcasecmp(sc1, interface_commands[ i ]) == 0)
          command_number = i;
    }
    else
	command_number = 14;

    switch (command_number)
    {
      case 0:
        exit = 1;
        break;
      case 1:
        center_img(&st[ 0 ], user.iNX);
        break;
      case 2:
        display_oifits_info(&data);
        break;
      case 3:
        goto RESTART;
        break;
      case 4:
        dispuv(dispdev);
        break;
      case 5:
        display_command_help();
        break;
      case 6:
        if (disp_id == 0)
        {
          open_redisplay(&user, dispdev);
          redisplay(&st[ 0 ], &user, reconst_parameters.displaypower, contour);
          get_default_model(startmod, &st[ 0 ], user.iNX, 0.0, user.iNX * user.xyint, 0.0, user.iNY * user.xyint, user.xyint, disp_id);
          redisplay(startmod, &user, reconst_parameters.displaypower, contour);
          close_redisplay(&user);
          disp_id = 0;
        }
        else
        {
          open_redisplay(&user, dispdev);
          redisplay(&st[ 0 ], &user,  reconst_parameters.displaypower, contour);
          get_default_model(startmod, &st[ 0 ], user.iNX, 0.0, user.iNX * user.xyint, 0.0, user.iNY * user.xyint, user.xyint, disp_id);
          redisplay(startmod, &user,  reconst_parameters.displaypower, contour);
        }
        goto RESTART;
        break;

      case 7:
        printf("Creating stddev map...\n");
        ncorr = 1;
        samples = 500;
        maps = malloc(samples * user.iNX * user.iNY * sizeof(float));
        errmap = malloc(user.iNX * user.iNY * sizeof(float));
        meanmap = malloc(user.iNX * user.iNY * sizeof(float));
        // Backup current image in st
        for (j = 0; j < user.iNX * user.iNY; j++)
        {
          errmap[ j ] = st[ j ];
        }
        for (i = 0; i < samples; i++)
        {
          MOVIE4(st, &ncorr, &kstat, &ntrnsx);
          fprintf(stderr, "Sample %d NTRNSX === %d kstat === %d \n", i + 1, ntrnsx, kstat);
          /* Store the current image */
          for (j = 0; j < user.iNX * user.iNY; j++)
          {
            maps[ j + i * user.iNX * user.iNY ] = st[ j ];
          }
        }
        /* Now compute the standard deviation map */

        for (j = 0; j < user.iNX * user.iNY; j++)
        {
          /* Mean */
          meanmap[ j ] = 0.;
          for (i = 0; i < samples; i++)
          {
            meanmap[ j ] += maps[ j + i * user.iNX * user.iNY ];
          }
          meanmap[ j ] /= ((float) samples);

          sum = 0.;
          for (i = 0; i < samples; i++)
          {
            sum += square(maps[ j + i * user.iNX * user.iNY ] - meanmap[ j ]);
          }
          //ERRMAP=image/standard deviation
          errmap[ j ] = sqrt(1. / ((float) samples - 1.) * sum);
        }

        ///* Reassign it to st ( for writefits purpose ) */
        for (j = 0; j < user.iNX * user.iNY; j++)
          st[ j ] = errmap[ j ];

        /* Display the final error map */
        redisplay(&errmap[ 0 ], &user,  reconst_parameters.displaypower, contour);
        free(maps);
        free(meanmap);
        free(errmap);
        break;

      case 8:
        printf("Creating average map...\n");
        ncorr = 1;
        samples = 50;
        meanmap = malloc(user.iNX * user.iNY * sizeof(float));
        for (j = 0; j < user.iNX * user.iNY; j++)
          meanmap[ j ] = 0.;
        for (i = 0; i < samples; i++)
        {
          MOVIE4(st, &ncorr, &kstat, &ntrnsx);
          fprintf(stderr, "Sample %d NTRNSX === %d kstat === %d \n", i + 1, ntrnsx, kstat);
          for (j = 0; j < user.iNX * user.iNY; j++)
            meanmap[ j ] += st[ j ];

        }

        ///* Reassign it to st ( for writefits purpose ) */
        for (j = 0; j < user.iNX * user.iNY; j++)
          st[ j ] = meanmap[ j ] / ((float) samples);

        /* Display the final error map */
        redisplay(&st[ 0 ], &user, reconst_parameters.displaypower, contour);
        free(meanmap);
        break;

      case 9:
        // filename = '!'+filename needed for CFITSIO to be able to overwrite if existing
        strcpy(reconst_parameters.fitsimage, sc2);
        if (write_fits_image(&st[ 0 ], &reconst_parameters, &status))
          printf("Error creating fits image %s.\n", sc2);
        break;

      case 10:
        sscanf(sc2, "%f", &reconst_parameters.displaypower);
        break;

      case 11:
        read_model(sc2, startmod, user.iNX, reconst_parameters.modelflux);
        redisplay(startmod, &user,reconst_parameters.displaypower, contour);
        goto RESTART;
        break;

      case 12:
        if (strcasecmp(sc2, "ON") == 0)
          reconst_parameters.verbose = 1;
        if (strcasecmp(sc2, "OFF") == 0)
          reconst_parameters.verbose = 0;
        break;

      case 13:
        if (strcasecmp(sc2, "ON") == 0)
        {
          open_redisplay(&user, dispdev);
          redisplay(&st[ 0 ], &user, reconst_parameters.displaypower, contour);
          disp_id = 1;
        }
        if (strcasecmp(sc2, "OFF") == 0)
        {
          close_redisplay(&user);
          disp_id = 0;
        }
        break;

      case 14:

        if (reconst_parameters.noui == 0)
        {
          sscanf(sc2, "%d", &si2);
          if (si2 < 0)
            si2 = 200000;
        }
        else
          si2 = 200000;

        stop = 0;
        i = 0;
        do
        {

          MEM4(st, &memrun, &entropy, &test, &chisq, &scale, &plow, &phigh, &pdev, &glow, &ghigh, &gdev, &omega, &alpha, &istat, &ntrans);

          total_iter++;
          i++;
          if (reconst_parameters.verbose == 1)
          {
            printf("Iteration %d Ntrans === %d istat === %d%d%d%d%d%d%d \n", total_iter, ntrans, istat % 2, istat / 2 % 2, istat / 4 % 2,
                istat / 8 % 2, istat / 16 % 2, istat / 32 % 2, istat / 64 % 2);
            printf("Entropy === %f  Chisq === %f Flux === %f Alpha === %f Omega === %f \n", entropy, chisq / (float) ndata, user.imflux, alpha, omega);
            printf("Logprob === %f Test === %f Good Measurements === %f Scale === %f \n", (plow + phigh) / 2., test, (glow + ghigh) / 2.,
                scale);
          }

          if (isnan(omega) > 0)
          {
            printf("Omega has been affected a NaN value.\n");
            printf("Model is probably inadequate.\n");
            break;
          }

          if (disp_id > 0)
            redisplay(&st[ 0 ], &user, reconst_parameters.displaypower, contour);

          if (((istat % 2) == 0) && ((istat / 2 % 2) == 0))
          {
            stop++;
            // Experimental -- recomputes the errors using the current visibilities once convergence is achieved
            //reset_memsys_dataspace( st, kb, &data, &reconst_parameters );
          }
          else
            stop = 0;
          if (i >= si2)
            stop = 6;
          if ((total_iter >= reconst_parameters.maxiter) && (reconst_parameters.maxiter != -1))
          {
            printf("Maximum iterations has been reached.\n");
            stop = 6;
          }

        } while (stop < 5);

        if (reconst_parameters.noui == 1)
        {
          if (write_fits_image(&st[ 0 ], &reconst_parameters, &status))
            printf("Error creating fits image %s.\n", sc2);
          exit = 1;
        }
        break;

      case 15: // Contours
        if (contour == 1)
          contour = 0;
        else
          contour = 1;
        break;

      case 16:

              break;

      default:
        printf("Unrecognized command\n");
        break;

    }
  } while (exit == 0);

  nfft_finalize(&user.p);
  free(st);
  free(startmod);
  free(user.current_cvis);
  free(user.current_derivcvis);
  EXIT: printf("\n\n");

  return SUCCESS;
}

int read_model( char* filename , float* mod , int N , float normalisation)
{
  int i, j;
  int nr;
  float ftmp;
  float flux = 0.0;
  char tmpch[ 20 ];
  float *tmpmod;
  int status = 0;
  int warnneg = 0;
  float minvalue = 1e-9;
  /* Create a temporary large model */
  tmpmod = malloc(N * N * sizeof(float));
  read_fits_image(filename, tmpmod, &nr, &ftmp, tmpch, tmpch, &status);

  if (status != 0)
  {
    printf("Error opening file %s\n", filename);
    return 0;
  }
  if (nr != N)
  {
    printf("Model image must have dimensions %d x %d pizels.\n", N, N);
    return 0;
  }


  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
        mod[ i * N + N - j - 1 ] = tmpmod[ i * N + j ];
	//	mod [i * N + j ]= tmpmod[i * N +j];
    }
  }
  free(tmpmod);

  for (i = 0; i < (N * N); i++)
   {
     if (mod[ i ] <= 0.0)
       {
	 warnneg = 1;
	 mod[ i ] = minvalue; /* set then to very small */
       }
     flux += mod[ i ];
   }

  for (i = 0; i < (N * N); i++)
    mod[ i ] *= normalisation / flux;

  if (warnneg == 1)
  {
    printf("Values <= 0 in model are not permitted. \n");
    printf("All have been replaced by %f.\n", minvalue);
  }

  return 1;
}

void uvchuck( float *st , int *kb , RECONST_PARAMETERS *reconst_parameters )
{
  int i;
  float uscale, vscale, ucut, vcut, ovsamp;
  float ftmp;

  /* Find maximum UV coord */
  float max_uv = 0.0;
  float min_uv = 1e20;
  for (i = 1; i < nuv; i++)
  {
      ftmp = u[i]*u[i] + v[i]*v[i];
      if (ftmp > max_uv)
        max_uv = ftmp;
      if (ftmp < min_uv)
        min_uv = ftmp;
  }
  max_uv = sqrt(max_uv);
  min_uv = sqrt(min_uv);
  printf("UV range:\t\t%8.0f - %8.0f wavelengths \n", min_uv, max_uv);

  /* Automatic choice of pixellation if xyint <= 0. UV cut off such that
   * the oversampling >= defined value ( OVERSAMPLING in constants.h ).
   * xyint in milliarcsec. uv in wavelengths.
   */
  printf("Array resolution:\t %f mas\n", 1. / (2.0 * MAS * max_uv));

  if (user.xyint <= 0.0)
  {
    user.xyint = 1.0 / (2.0 * MAS * OVERSAMPLING * max_uv);
    printf("Pixel size:\t\tAutomatic, %f mas\n", user.xyint);
    reconst_parameters->xyint = user.xyint;
    printf("Recommended size:\t %d pixels\n", (int) (powf(2.0, ceil(log(max_uv / min_uv * 2. * OVERSAMPLING + 1.) / log(2.) ) )));
  }
  else
    printf("Pixel size:\t\tUser defined, %f mas\n", user.xyint);

  ovsamp = 1.0 / (max_uv * user.xyint * MAS);

  printf("Image width:\t\t%d pixels, %f mas\n", user.iNX, (float) user.iNX * user.xyint);
  printf("Pix/fastest fringe:\t%f\n", ovsamp);


}

int set_memsys_dataspace( float *st , int *kb  , RECONST_PARAMETERS *reconst_parameters )
{
  // Note: errors bar in the data should respect the strict OIFITS definition
  // ( cf PASP paper ) - in particular, closures may be >180
  int i;
  int warn_extrapolation = 0;
  float err_pow, err_rad, err_tan, err_abs, err_phi=0;
  float pow1, powerr1, pow2, powerr2, pow3, powerr3, sqamp1, sqamp2, sqamp3, sqamperr1, sqamperr2, sqamperr3;
  printf("Loading data into memory\n");

  for (i = 0; i < nv2; i++)
  {
    st[ kb[ 20 ] + i ] = data[ i ];
    st[ kb[ 21 ] + i ] = 1.0 / data_err[ i ] ;
 }

  for (i = 0; i < nt3; i++)
  {

    // if  ( (reconst_parameters->forced_extrapolation == 1) || ((data_err[nv2 + i ] <= 0.) || (data_err[nv2 + i ] > (infinity - 1)))) // Missing triple amplitudes
    //{
      // check all pows exist  if ((t3in1[ i ] < nv2) && (t3in2[ i ] < nv2) && (t3in3[ i ] < nv2))

      // if corresponding powerspectrum points are available
    //      {
	// Derive pseudo-triple amplitudes from powerspectrum data
        // First select the relevant powerspectra
    //  pow1 = data[ t3in1[ i ] ];
    //   powerr1 = data_err[[ t3in1[ i ] ];
    //   pow2 = data[ t3in2[ i ] ];
    //   powerr2 = data_err[[ t3in2[ i ] ];
    //   pow3 = data[ t3in3[ i ] ];
    //   powerr3 = data_err[[ t3in3[ i ] ];
        // Derive unbiased visibility amplitudes + noise variance
    //        sqamp1 = (pow1 + sqrt(square(pow1) + 2.0 * square(powerr1))) / 2.;
    //   sqamperr1 = 1. / (1. / sqamp1 + 2. * (3. * sqamp1 - pow1) / square(powerr1));
    //   sqamp2 = (pow2 + sqrt(square(pow2) + 2.0 * square(powerr2))) / 2.;
    //  sqamperr2 = 1. / (1. / sqamp2 + 2. * (3. * sqamp2 - pow2) / square(powerr2));
    //  sqamp3 = (pow3 + sqrt(square(pow3) + 2.0 * square(powerr3))) / 2.;
    //  sqamperr3 = 1. / (1. / sqamp3 + 2. * (3. * sqamp3 - pow3) / square(powerr3));

        // And form the triple amplitude statistics
    //  data[nv2 + i ] = sqrt(sqamp1 * sqamp2 * sqamp3);
    //  data_err[nv2 + i ] = fabs(data[nv2 + i ]) * sqrt(sqamperr1 / sqamp1 + sqamperr2 / sqamp2 + sqamperr3 / sqamp3);
	//   printf("Triple amplitude %d extrapolated from powerspectrum data T = %f \t E_T = %f: \r", i, data[nv2 + i ], data_err[nv2 + i ]);
    //	if(warn_extrapolation == 0)
    //	  {
    //	    warn_extrapolation = 1;
    //	    printf("WARNING: at least one triple amplitude is missing and has been extrapolated from powerspectrum values\n");
    //	  }
    //  }

    //  else // missing powerspectrum points -> cannot extrapolate bispectrum
    //  {
    //    printf("WARNING: triple amplitude extrapolation from powerspectrum failed because of missing powerspectrum\n");
    //	data[nv2 + i ] = 1.0;
    //    data_err[nv2 + i ] = infinity;
    //  }
    // }

    if(reconst_parameters->t3errtype == 1 )
      {
	// Full elliptic approximation - Initially based on Meimon 2009 appendix E 2-3
	// See also WISARD manual

	if (i == 0)
	  printf("Bispectrum noise:\tImproved elliptic approximation \n");


	err_abs = data_err[nv2 + i ];
	err_phi = data_err[nv2 + nt3amp + i ];
	err_rad =  sqrt(0.5 * square( err_abs )   * square(1.+ exp(-2.*square(err_phi)))
		      + 0.5 * square( data[nv2 + i ]) * square(1.- exp(-square(err_phi)))  );

	err_tan =  sqrt(0.5 * square( err_abs )   * square(1.- exp(-2.*square(err_phi)))
		      + 0.5 * square( data[nv2 + i ]) * (1.- exp(-2.*square(err_phi))) );

      }
    else
      {
	// Approximation - 1st order
	if (i == 0)
	  printf("Bispectrum noise:\tClassic elliptic approximation \n");
	err_rad = data_err[nv2 + i ] ;
	err_tan = fabs(data[nv2 + i ] * data_err[nv2 + nt3amp + i ] * M_PI / 180.0);
      }

    //
    // Set bispectrum errors
    //
    st[ kb[ 21 ] + nv2 + 2 * i ] = 1.0 / err_rad;
    st[ kb[ 21 ] + nv2 + 1 + 2 * i ] = 1.0 / err_tan;

    // Set Bispectrum data, rotated by exp( -i*closure_data )
    //
    if(reconst_parameters->t3errtype == 1 )
	// Improved elliptic approximation
      st[ kb[ 20 ]+ nv2 + 2 * i ] = data[nv2 + i ] * (2.- exp( - 0.5*square(err_phi) ))  ;
   else
	// Approximation 1st order
        st[ kb[ 20 ] + nv2 + 2 * i ] = data[nv2 + i ];

    st[ kb[ 20 ] + nv2 + 1 + 2 * i ] = 0.0;

  }
  printf("Data loaded in memory\n");
  return SUCCESS;
}

// USAVE and UREST are needed for linking with MEMSYS library

void __stdcall USAVE( int *ints , int *nints , float *reals , int *nreals )
{// Routine for saving crucial MemSys areas and scalars on disk
}

void __stdcall UREST( int *ints , int *nints , float *reals , int *nreals )
{// Routine for restoring crucial MemSys areas and scalars from disk
}

void display_oifits_info( )
{
  float ratio, nmin1 = 0.0, nmin2 = 0.0, nmin3 = 0.0, nmax1 = 1e10, nmax2 = 1e10, nmax3 = 1e10, tot1 = 0., tot2 = 0., tot3 = 0.;
  char choice[ 5 ], tempstr;
  int ii;
  printf("Display powerspectrum data ? y/[ N ] ");
  if(fgets(choice, sizeof(choice), stdin) ==NULL)
    printf("Error getting answer\n");
  tempstr = choice[ 0 ];
  for (ii = 0; ii < nv2; ii++)
  {
    ratio = fabs(data[ ii ] / data_err[ ii ]);
    if ((tempstr == 'y') || (tempstr == 'Y'))
      printf("N: %d Amp: %9f +/-%9f | S/N %3.1f\n", ii, data[ ii ], data_err[ ii ], ratio);
    if (ii == 1)
    {
      nmax1 = ratio;
      nmin1 = ratio;
    }
    if(ii > 1)
      {
	if (ratio > nmax1)
	  nmax1 = ratio;
	if (ratio < nmin1)
	  nmin1 = ratio;
	tot1 += ratio / nv2;
      }
  }
  printf("Display bispectrum data ? y/[ N ] ");
  if(fgets(choice, sizeof(choice), stdin) ==NULL)
    printf("Error getting answer\n");

  tempstr = choice[ 0 ];
  for (ii = 0; ii < nt3; ii++)
  {
    ratio = fabs(data[nv2 + ii ] / data_err[nv2 + ii ]);
    if ((tempstr == 'y') || (tempstr == 'Y'))
      printf("N: %d Amp : %9f +/-%9f | S/N %3.1f | Phs : %9f +/-%9f\n", ii, data[nv2 + ii ], data_err[nv2 + ii ], ratio,
          data[nv2 + nt3amp + ii ], data_err[nv2 + nt3amp + ii ]);
    if (ii == 0)
    {
      nmin2 = ratio;
      nmax2 = ratio;
      nmin3 = data_err[nv2 + nt3amp + ii ];
      nmax3 = data_err[nv2 + nt3amp + ii ];
    }
    if (ratio < nmin2)
      nmin2 = ratio;
    if (ratio > nmax2)
      nmax2 = ratio;
    if (data_err[nv2 + nt3amp + ii ] < nmin3)
      nmin3 = data_err[nv2 + nt3amp + ii ];
    if (data_err[nv2 + nt3amp + ii ] > nmax3)
      nmax3 = data_err[nv2 + nt3amp + ii ];
    tot2 += ratio / nt3;
    tot3 += data_err[nv2 + nt3amp + ii ] / nt3;
  }

  printf("SNR |   v2    |  t3amp    | t3phs    |\n");
  printf("Min | %8.1f | %8.1f | %8.3f |\n", nmin1, nmin2, nmin3);
  printf("Avg | %8.1f | %8.1f | %8.3f |\n", tot1, tot2, tot3);
  printf("Max | %8.1f | %8.1f | %8.3f |\n", nmax1, nmax2, nmax3);

}

void set_model( RECONST_PARAMETERS* reconst_parameters , float *startmod )
{

  if (strcmp(reconst_parameters->priorfile, "") != 0)
  {
    if (reconst_parameters->verbose == 1)
      printf("Prior read from file = %s\n", reconst_parameters->priorfile);
    read_model(reconst_parameters->priorfile, startmod, user.iNX, reconst_parameters->modelflux);
  }
  else
  {
  float flux = 0.0;
    int i, j;
    switch (reconst_parameters->modeltype)
    {
      case 0:
      {
        // Flat prior - generally does not work as the centering of the reconstruction depends on the prior !
        if (reconst_parameters->verbose == 1)
          printf("Flat prior, total flux: %8.3f\n", reconst_parameters->modelflux);
        for (i = 0; i < user.iNX; i++)
          for (j = 0; j < user.iNY; j++)
            startmod[ j * user.iNX + i ] = reconst_parameters->modelflux / ((float) user.iNY * (float) user.iNX);

        break;
      }

      case 1:
      {
        // Centered Dirac
        if (reconst_parameters->verbose == 1)
          printf("Dirac, flux: %8.3f\n", reconst_parameters->modelflux);

        for (i = 0; i < user.iNX; i++)
          for (j = 0; j < user.iNY; j++)
            startmod[ j * user.iNX + i ] = 1e-8;

        startmod[ (user.iNY * (user.iNX + 1)) / 2 ] = reconst_parameters->modelflux;

        break;
      }

      case 2:
      {
        // Centered Uniform disk
        if (reconst_parameters->verbose == 1)
          printf("Uniform disk, Radius:%f mas, flux:%f\n", reconst_parameters->modelwidth, reconst_parameters->modelflux);
        flux = 0.0;
        float rsq;
        for (i = 0; i < user.iNX; i++)
        {
          for (j = 0; j < user.iNY; j++)
          {
            rsq = square((float) (user.iNX / 2 - i)) + square((float) (user.iNY / 2 - j));

            if (sqrt(rsq) <= (reconst_parameters->modelwidth / user.xyint))
              startmod[ j * user.iNX + i ] = 1.;
            else
              startmod[ j * user.iNX + i ] = 1e-8;
            flux += startmod[ j * user.iNX + i ];
          }
        }

        for (i = 0; i < user.iNX * user.iNY; i++)
          startmod[ i ] *= reconst_parameters->modelflux / flux;

        break;
      }

      case 3:
      {
        // Centered Gaussian
        float sigma = reconst_parameters->modelwidth / (2. * sqrt(2. * log(2.)));
        if (reconst_parameters->verbose == 1)
          printf("Gaussian, FWHM:%f mas, sigma:%f mas, flux:%f\n", reconst_parameters->modelwidth, sigma, reconst_parameters->modelflux);
        flux = 0.0;
        for (i = 0; i < user.iNX; i++)
        {
          for (j = 0; j < user.iNY; j++)
          {
            startmod[ j * user.iNX + i ] = exp(-(square((float) user.iNX / 2 - i) + square((float) user.iNY / 2 - j)) / (2. * square(sigma
                / user.xyint)));
            // Fix problem with support at zero
            //if ( startmod[ j*user.iNX+i ]<1e-8 ) startmod[ j*user.iNX+i ]=1e-8;
            flux += startmod[ j * user.iNX + i ];
          }
        }
        for (i = 0; i < user.iNX * user.iNY; i++)
          startmod[ i ] *= reconst_parameters->modelflux / flux;
        break;
      }

      case 4:
      {
        //Centered Lorentzian
        float sigma = reconst_parameters->modelwidth;
        if (reconst_parameters->verbose == 1)
          printf("Lorentzian, FWHM:%f mas, sigma:%f mas, flux:%f\n", reconst_parameters->modelwidth, sigma, reconst_parameters->modelflux);
        flux = 0.0;
        for (i = 0; i < user.iNX; i++)
        {
          for (j = 0; j < user.iNY; j++)
          {
            startmod[ j * user.iNX + i ] = sigma / (square(sigma) + square((float) user.iNX / 2 - i) + square((float) user.iNY / 2 - j));
            flux += startmod[ j * user.iNX + i ];
          }
        }
        for (i = 0; i < user.iNX * user.iNY; i++)
          startmod[ i ] *= reconst_parameters->modelflux / flux;

        break;
      }

    }

  }
}

float square( float n )
{
  return n * n;
}

int write_fits_image( float* img , RECONST_PARAMETERS* reconst_parameters , int* status )
{
  fitsfile *fptr;
  int i, j;
  long fpixel = 1, naxis = 2, nelements;
  long naxes[ 2 ];
  float* flipimg;
  char fitsimage[ 100 ];
  for (i = 0; i < 100; i++)
    fitsimage[ i ] = '\0';
  /*Initialise storage*/
  naxes[ 0 ] = (long) reconst_parameters->xydim;
  naxes[ 1 ] = (long) reconst_parameters->xydim;
  nelements = naxes[ 0 ] * naxes[ 1 ];

  if (strcmp(reconst_parameters->fitsimage, "") != 0)
  {
    strcpy(fitsimage, "!");
    strcat(fitsimage, reconst_parameters->fitsimage);
  }
  else
    strcpy(fitsimage, "!output.fits");

  // Flip the image so that East = Left
    flipimg = malloc(reconst_parameters->xydim * reconst_parameters->xydim * sizeof(float));
    for (i = 0; i < reconst_parameters->xydim; i++)
  for (j = 0; j < reconst_parameters->xydim; j++)
    flipimg[ i * reconst_parameters->xydim + reconst_parameters->xydim - j - 1 ] = img[ i * reconst_parameters->xydim + j ];

  /*Create new file*/
  if (*status == 0)
    fits_create_file(&fptr, fitsimage, status);

  /*Create primary array image*/
  if (*status == 0)
    fits_create_img(fptr, FLOAT_IMG, naxis, naxes, status);
  /*Write a keywords (datafile, target, image pixelation) */
  if (*status == 0)
    fits_update_key(fptr, TSTRING, "DATAFILE", reconst_parameters->datafile, "Data File Name", status);
  if (*status == 0)
    fits_update_key(fptr, TSTRING, "TARGET", reconst_parameters->targetname, "Target Name", status);
  if (*status == 0)
    fits_update_key(fptr, TFLOAT, "PIXSIZE", &reconst_parameters->xyint, "Pixelation (mas)", status);
  if (*status == 0)
    fits_update_key(fptr, TINT, "WIDTH", &reconst_parameters->xydim, "Size (pixels)", status);
   if(*status == 0)
     fits_update_key(fptr, TFLOAT, "LOGZLOW", &plow,"Evidence (low)", status);
   if(*status == 0)
     fits_update_key(fptr, TFLOAT, "LOGZHI", &phigh,"Evidence (high)", status);
   if(*status == 0)
     fits_update_key(fptr, TFLOAT, "ALPHA", &alpha, "Regularization hyperparameter", status);
   if(*status == 0)
     fits_update_key(fptr, TINT, "NITER", &total_iter,"Number of iterations", status);

  /*
   if(*status == 0)fits_update_key(fptr, TFLOAT, "CHI2", &chisq, "Chi2", status);
   float logprob = (plow+phigh)/2.;


   if(*status == 0)fits_update_key(fptr, TFLOAT, "OMEGA", &omega, "Reconstruction success indicator", status);

   */


  /*Write image*/
  if (*status == 0)
    fits_write_img(fptr, TFLOAT, fpixel, nelements, &flipimg[ 0 ], status);

  /*Close file*/
  if (*status == 0)
    fits_close_file(fptr, status);

  /*Report any errors*/
  fits_report_error(stderr, *status);

  /*Error handling*/
  free(flipimg);
  return *status;
}

void __stdcall VMEMEX( float *image , float *model_data )
{
  // Visible-to-Data transform
  // Takes present model image and calculates the appropriate non-linear
  // mock data values. (Powerpectrum and bispectrum points)
  float complex vtemp, V0ab, V0bc,V0ca;
  float complex* cvis = user.current_cvis;
  
  
  user.imflux = 0.0;  
  for (int i = 0; i < user.iNX * user.iNX; i++)
    {
      user.p.f_hat[ i ] = image[ i ] + I * 0.0;
      user.imflux += image[ i ];
    }
  nfft_trafo(&user.p);
  
  for (int uu = 0; uu < nuv; uu++)
    user.current_cvis[ uu ] = user.p.f[ uu ];

  // Powerpesctrum
  for (int i = 0; i < nv2; i++)
      model_data[ i ] = abs2(user.current_cvis[ i ]);

  //printf("Check - powerspectrum 0 == %f\n", model_data[0]);
  // Bispectrum
  for (int i = 0; i < nt3; i++)
  {
      V0ab = cvis[ t3in1[ i ] ];
      V0bc = cvis[ t3in2[ i ] ];
      V0ca = cvis[ t3in3[ i ] ];
      vtemp = V0ab * V0bc * conj(V0ca) * user.data_phasor[ i ];

      model_data[ nv2 + 2 * i ] = creal(vtemp);
      model_data[ nv2 + 2 * i + 1 ] = cimag(vtemp);

  }

}

void __stdcall VOPUS( float *dh , float *dd )
{
  // Visible-to-Data differential transform
  // Note current_derivcvis is dVisibility, dh is dImage, dd is dData
  float complex vtemp;
  float complex V0ab, V0bc, V0ca;
  float complex VISab, VISbc, VISca;
  float complex* cvis = user.current_cvis;
  float complex* dcvis = user.current_derivcvis;

  
  for (int i = 0; i < user.iNX * user.iNX; i++)
    user.p.f_hat[ i ] = dh[ i ] + I * 0.0;

  nfft_trafo(&user.p);

  for (int u = 0; u < nuv; u++)
        user.current_derivcvis[ u ] = user.p.f[ u ];

  /* Powerspectrum */
  for (int i = 0; i < nv2; i++)
      dd[ i ] = 2.0 *
	(creal(dcvis[ i ]) * creal(cvis[ i ]) + cimag(dcvis[ i ]) * cimag(cvis[ i ]));


  /* Bispectrum */
  for (int i = 0; i < nt3; i++)
  {
      V0ab = cvis[ t3in1[ i ] ];
      V0bc = cvis[ t3in2[ i ] ];
      V0ca = cvis[ t3in3[ i ] ];
      VISab = dcvis[t3in1[ i ]  ];
      VISbc = dcvis[ t3in2[ i ] ];
      VISca = dcvis[ t3in3[ i ] ];
      /* differential response calculation */
      vtemp = user.data_phasor[ i ] * (VISab * V0bc * conj(V0ca) + VISbc * conj(V0ca) * V0ab + conj(VISca) * V0ab * V0bc);
      dd[ nv2 + 2 * i ] = creal(vtemp);
      dd[ nv2 + 2 * i + 1 ] = cimag(vtemp);

  }

}

void __stdcall VTROP( float *dh , float *dd )
{ // Data-to-Visible differential transform

  float complex V0ab,V0bc,V0ca;
  float complex VISab,VISbc,VISca;
  float complex t3;
  float complex* cvis = user.current_cvis;
  float complex* dcvis = user.current_derivcvis;

  for (int i = 0; i < nuv; i++)
    dcvis[ i ] = 0.0;

  /* Powerspectrum contribution */
  for (int i = 0; i < nv2; i++)
      dcvis[ i ] += cvis[ i ] * 2.0 * dd[ i ];

  /* Bispectrum contribution */
  for (int i = 0; i < nt3; i++)
  {
      /* previous visibilities */
      V0ab = conj(cvis[t3in1[i]]);
      V0bc = conj(cvis[t3in2[i]]);
      V0ca = cvis[t3in3[i]];

      /* input bispectrum differential */
      t3 = (dd[ nv2 + 2 * i ] + I * dd[ nv2 + 2 * i + 1 ]) * conj(user.data_phasor[ i ]);

      /* differential response calculation */
      VISab = t3 * V0bc * V0ca;
      VISbc = t3 * V0ca * V0ab;
      VISca = conj(t3 * V0ab * V0bc);

      /* Add to differential response */
      dcvis[ t3in1[i] ] += VISab;
      dcvis[ t3in2[i] ] += VISbc;
      dcvis[ t3in3[i] ] += VISca;

  }

  for (int uu = 0; uu < nuv; uu++)
    user.p.f[ uu ] = dcvis[ uu ];
  
  nfft_adjoint(&user.p);
  
  for (int ii = 0; ii < user.iNX * user.iNX; ii++)
      dh[ ii ] = creal(user.p.f_hat[ ii ]);

}


//void op_trop_check( int npix , int ndata )
//{
//  int i, jj;
//  float *hid1, *vis1, *hid2, *vis2, *temp;
//  float sumh, sumv;
//  float precision;
//  unsigned Rand[ 113 ];
//  RanInit(Rand, 1);// -1 for time based, +1 for non time based
//  printf("Testing Opus/tropus consistency with %d pixels and %d data\n", npix, ndata);
//
//  hid1 = malloc(npix * sizeof(float));
//  vis1 = malloc(ndata * sizeof(float));
//  hid2 = malloc(npix * sizeof(float));
//  vis2 = malloc(ndata * sizeof(float));
//
//  temp = malloc(ndata * sizeof(float));
//
//  for (jj = 0; jj < 10; jj++)
//  {
//
//    for (i = 0; i < npix; i++)
//    {
//      hid1[ i ] = Ranfloat(Rand);
//      hid2[ i ] = 0.0;
//    }
//    for (i = 0; i < ndata; i++)
//    {
//      vis1[ i ] = Ranfloat(Rand);
//      vis2[ i ] = 0.0;
//    }
//
//    VMEMEX(hid1, temp);
//    VOPUS(hid1, vis2);
//    VTROP(hid2, vis1);
//    sumh = 0.0;
//    sumv = 0.0;
//    for (i = 0; i < npix; i++)
//      sumh += hid1[ i ] * hid2[ i ];
//    for (i = 0; i < ndata; i++)
//      sumv += vis1[ i ] * vis2[ i ];
//    precision = fabs((sumh - sumv) / sumh);
//    printf("OPTROP Test %d Precision : %e ( SUMH : %f SUMV : %f ) \n", jj, precision, sumh, sumv);
//  }
// free(vis1);
// free(vis2);
// free(hid1);
// free(hid2);
// free(temp);
//
//
//}

float abs2(float complex cp)
{
  return creal(cp)*creal(cp)+cimag(cp)*cimag(cp);
}

void center_img( float* img , int N ) // recenter an image (avoid this by choosing a centered starting prior)
{
  int i, k, l;
  float mx = 0.0, my = 0.0, min = 1e9, max = 0.0;
  float flux = 0.0;
  float* tmp;
  int cx, cy;
  int dx, dy;

  tmp = malloc(N * sizeof(float));

  for (i = 0; i < N; i++)
  {
    for (k = 0; k < N; k++)
    {
      if (max < img[ i * N + k ])
        max = img[ i * N + k ];
      if (min > img[ i * N + k ])
        min = img[ i * N + k ];
    }
  }

  /* Find position of CofM on thresholded image */
  for (i = 0; i < N; i++)
  {
    for (k = 0; k < N; k++)
    {
      if (img[ i * N + k ] > max / 10.)
      {
        mx += img[ i * N + k ] * (float) i;
        my += img[ i * N + k ] * (float) k;
        flux += img[ i * N + k ];
      }
    }
  }

  cx = (int) (mx / flux + 0.5);
  cy = (int) (my / flux + 0.5);

  dx = N / 2 - cx;
  dy = N / 2 - cy;

  printf("Moving image by ( %d, %d ) pixels.\n", dx, dy);

  /* X */
  if (dx > 0)
  {
    for (l = 0; l < dx; l++)
    {
      for (k = 0; k < N; k++)
        tmp[ k ] = img[ (N - 1) * N + k ];

      for (i = 0; i < (N - 1); i++)
      {
        for (k = 0; k < N; k++)
        {
          img[ (N - 1 - i) * N + k ] = img[ (N - 2 - i) * N + k ];
        }
      }

      for (k = 0; k < N; k++)
        img[ k ] = tmp[ k ];
    }
  }
  if (dx < 0)
  {
    for (l = 0; l < (-dx); l++)
    {
      for (k = 0; k < N; k++)
        tmp[ k ] = img[ k ];

      for (i = 0; i < (N - 1); i++)
      {
        for (k = 0; k < N; k++)
        {
          img[ i * N + k ] = img[ (i + 1) * N + k ];
        }
      }

      for (k = 0; k < N; k++)
        img[ (N - 1) * N + k ] = tmp[ k ];

    }
  }

  /* Y */
  if (dy > 0)
  {
    for (l = 0; l < dy; l++)
    {
      for (k = 0; k < N; k++)
        tmp[ k ] = img[ k * N + (N - 1) ];

      for (i = 0; i < (N - 1); i++)
      {
        for (k = 0; k < N; k++)
        {
          img[ i * N + (N - 1 - k) ] = img[ i * N + (N - 2 - k) ];
        }
      }

      for (k = 0; k < N; k++)
        img[ k * N ] = tmp[ k ];
    }
  }
  if (dy < 0)
  {
    for (l = 0; l < (-dy); l++)
    {
      for (k = 0; k < N; k++)
        tmp[ k ] = img[ k * N ];

      for (i = 0; i < (N - 1); i++)
      {
        for (k = 0; k < N; k++)
        {
          img[ i * N + k ] = img[ i * N + (k + 1) ];
        }
      }

      for (k = 0; k < N; k++)
        img[ k * N + (N - 1) ] = tmp[ k ];

    }
  }
  free(tmp);
}

int compare(char *A,char *B);
int commandline(RECONST_PARAMETERS* reconst_parameters, int argc, char** argv);
void reset_default_parameters(RECONST_PARAMETERS* reconst_parameters);
int get_input(RECONST_PARAMETERS* reconst_parameters, int argc, char** argv );
int interactivemode(RECONST_PARAMETERS* reconst_parameters);

// GET INPUT : switch between commandline and interactive mode
int get_input(RECONST_PARAMETERS* reconst_parameters, int argc, char** argv )
{
	int exit = 0;
	exit = commandline(reconst_parameters, argc, argv);
	return exit;
}

int commandline(RECONST_PARAMETERS* reconst_parameters, int argc, char** argv)
{
  printf("\n\n**********  BSMEM v2.0   ******************\n");
  int i;
  int help = 0;
  // DISPLAY HELP
  if ( (argc < 2) || (strcmp(argv[1],"-h") == 0) || (strcmp(argv[1],"-help") == 0)|| (strcmp(argv[1],"--help") == 0))
    help=1;
  
  // GET COMMANDLINE INPUT PARAMETERS
  reset_default_parameters(reconst_parameters);
  for(i=1 ; i < argc; i+=2)
    {
      if(strcmp(argv[i],"-d") == 0)
	sscanf(argv[i+1], "%s", reconst_parameters->datafile);
      else if (strcmp(argv[i],"-p") == 0)
	sscanf(argv[i+1],"%f", &reconst_parameters->xyint);
      else if (strcmp(argv[i],"-it") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->maxiter);
      else if (strcmp(argv[i],"-a") == 0)
	sscanf(argv[i+1],"%f", &reconst_parameters->alpha);
      else if (strcmp(argv[i],"-e") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->entropy);
      else if (strcmp(argv[i],"-r") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->regularization);
      else if (strcmp(argv[i],"-mf") == 0)
	sscanf(argv[i+1],"%f", &reconst_parameters->modelflux);
      else if (strcmp(argv[i],"-mw") == 0)
	sscanf(argv[i+1],"%f", &reconst_parameters->modelwidth);
      else if (strcmp(argv[i],"-mt") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->modeltype);
      else if (strcmp(argv[i],"-p") == 0)
	sscanf(argv[i+1],"%f", &reconst_parameters->xyint);
      else if (strcmp(argv[i],"-wavmin") == 0)
	  sscanf(argv[i+1],"%f", &reconst_parameters->minband);
      else if (strcmp(argv[i],"-wavmax") == 0)
	  sscanf(argv[i+1],"%f", &reconst_parameters->maxband);
      else if (strcmp(argv[i],"-w") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->xydim);
      else if (strcmp(argv[i],"-vb") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->verbose);
      else if(strcmp(argv[i],"-f") == 0)
	sscanf(argv[i+1], "%s", reconst_parameters->fitsimage);
      else if(strcmp(argv[i],"-sf") == 0)
	    sscanf(argv[i+1], "%s", reconst_parameters->priorfile);
      else if (strcmp(argv[i],"-berr") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->t3errtype);
      else if (strcmp(argv[i],"-ferr") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->fluxerr);
      else if(strcmp(argv[i],"-noui") == 0)
	reconst_parameters->noui = 1;
      else if(strcmp(argv[i],"-forcext") == 0)
	reconst_parameters->forced_extrapolation = 1;
      else if(strcmp(argv[i],"-ver") == 0)
	{
	  printf("***********************************************\n");
	  printf("*          B S M E M  v2.0                    *\n");
	  printf("*                                             *\n");
	  printf("* Main developper  : Fabien Baron             *\n");
	  printf("* OIFITS library   : John Young               *\n");
	  printf("*                                             *\n");
	  printf("***********************************************\n");
	}
      else if(strcmp(argv[i],"-h") == 0)
	help = 1;
      else
	{
	  printf("A parameter was not recognized: %s \n", argv[i]);
	  printf("Type bsmem -h to display the command line help.\n");
	  return  1;
	}
    }      
      
  if( help == 1)
    {
      printf("Usage: bsmem -d OIFITSfile [-f outputimagefile -p pixellation -w imagewidth ...] \n");
      printf("Example: './bsmem -d data.oifits -p 0.1 -w 128' \n\n");
      printf("-h:\t\t Display this information.\n");
      printf("-d:\t\t OIFITS file containing the visibility data.\n");
      printf("-f:\t\t FITS file to output the reconstructed image.\n");
      printf("-sf:\t\t Starting image or prior file. Overrides the -mt command.\n");
      printf("-mt:\t\t Model/prior image type. \n");
      printf("\t\t   0 : Flat prior.\n");
      printf("\t\t   1 : Dirac, centered in the FOV.\n");
      printf("\t\t   2 : Uniform disk.\n");
      printf("\t\t   3 : Gaussian.\n");
      printf("\t\t   4 : Lorentzian.\n");
      printf("-mw:\t\t Model witdth (Gaussian and Uniform Disk only). \n");
      printf("-mf:\t\t Total flux of the model. \n");
      printf("-p:\t\t Size of a pixel (in mas). Set to 0 for automatic.\n");
      printf("-w:\t\t Width (in pixels) of the reconstructed image.\n");
      printf("-e:\t\t Entropy functional.\n\t\t   1: Gull-Skilling entropy.\n\t\t   4: Quadratic regularization.\n");
      printf("-r:\t\t Regularization hyperparameter evaluation. \n");
      printf("\t\t   1: Classic Bayesian with known noise scale.\n\t\t   2: Classic Bayesian with simultaneous estimation of the noise scale.\n");
      printf("\t\t   3: Uses the user-provided ALPHA value.\n\t\t   4: Historic Chisquared=Ndata (fast, but non-Bayesian).\n");
      printf("-a:\t\t Regularization hyperparameter (-r 3 to enable) \n");
      printf("-it:\t\t Number of maximum iterations (-1 to disable)\n");
      printf("-berr:\t\t Type of errors on the rotated bispectrum vector\n");
      printf("\t\t   0 : Classic elliptic approximation, 1st order only (default)\n");
      printf("\t\t   1 : 'Full' elliptic approximation (Meimon et al.).\n");
      printf("-ferr:\t\t Error on the zero flux powerspectrum (default = 1e-4).\n");
      printf("-wavmin:\t\t Lower wavelength value for data selection purpose (in nm).\n");
      printf("-wavmax:\t\t Higher wavelength value for data selection purpose (in nm).\n");
      printf("-forcext:\t\t Force triple amplitude computation from powerspectrum data.\n");
      printf("-vb:\t\t Verbose on/off.\n");
      printf("-noui:\t\t No user interface, BSMEM exits immediately after the reconstruction.\n");
      printf("-ver:\t\t Version info.\n");
      return  1;
    }
  
  if(reconst_parameters->datafile[0]=='\0')
    {
      printf("No input datafile -> Exiting... \n");
      return 1;
    }
  
  return  0;
}

void display_command_help()
{
	printf( "\n AVAILABLE COMMANDS:\n" );
	printf( "\n EXIT:\t Exit BSMEM.\n" );
	printf( "\n Iteration commands \n\n" );
	printf( " DO n :\t n > 0: perform n iterations of MEMSYS\n");
	printf( "       \t n < 0: iterates until stopping criterion is met.\n" );
	printf( "\n Model commands \n\n" );
	printf( " STARTMOD:\t start a model.\n" );
	printf( " READMOD file:\t read a prior image from file.\n" );
	printf( " RESTART:\t restart iterations from the beginning.\n" );
	printf( "\n Graphics display commands \n\n" );
	printf( " SCALE x:\t set the powerscale x of the graphic display.\n" );
	printf( " CONTOUR:\t enable/disable the contour overplot.\n");
	printf( " REDISP ON/OFF:\t when ON, a window is updated with the current image estimate.\n" );
	printf( " CENTER:\t shift center of mass of the current image to the FOV center.\n" );
	printf( " ERRMAP:\t compute the SNR map.\n" );
	printf( "\n Data info \n\n" );
	printf( " SNR:\t\t display data information and SNR statistics.\n" );
	printf( " UV:\t\t display UV coverage of data set.\n" );
	printf( "\n Save commands \n\n" );
	printf( " WRITEFITS filename.fits:\t saves the current image as 'filename.fits'.\n");
	printf( "\n\n" );
}




void reset_default_parameters(RECONST_PARAMETERS* reconst_parameters)
{
	int i;
	for(i = 0; i < 100; i++) reconst_parameters->datafile[i]='\0';
	reconst_parameters->t3io = 0;
	reconst_parameters->v2io = 0;
	reconst_parameters->modelwidth = 10.0;
	reconst_parameters->modelflux = 0.01;
	reconst_parameters->modeltype = 3;
	reconst_parameters->blur = 0.0;
	reconst_parameters->xyint = -1.0;
	reconst_parameters->box = 0;
	reconst_parameters->xydim = 128;
	reconst_parameters->entropy = 1;
	reconst_parameters->regularization = 4;
	reconst_parameters->maxiter = 200;
	reconst_parameters->alpha = 100.;
	reconst_parameters->verbose = 1;
	reconst_parameters->noui = 0;
	reconst_parameters->t3errtype = 0;
	reconst_parameters->fluxerr = 1e-2;
	reconst_parameters->forced_extrapolation = 0;
	reconst_parameters->minband = -1. ;
	reconst_parameters->maxband = -1. ;

	for(i = 0; i < 100; i++) reconst_parameters->fitsimage[i]='\0';
	strcpy( reconst_parameters->fitsimage,  "output.fits");
	for(i = 0; i < 100; i++) reconst_parameters->priorfile[i]='\0';
}


int s_equals(char *A, char *B)
{
	/*
	* Extracts string from B to A.
	* B = "name"  --> A = name
	* B = name    --> A = name
	*/
	int i;

	if(B[0]=='\"')
	{
		i=1;
		do{
			A[i-1]=B[i];
			i++;
			if((B[i]=='\"')||(B[i]=='\0'))break;
		}while(1);
	}

	else
	{
		i = 0;
		do{
			A[i]=B[i];
			i++;
			if((B[i]=='\"')||(B[i]=='\0'))break;
		}while(1);
	}

	return 1;
}

int get_default_model(float* model, float* currentimage, int n, float xs, float xe, float ys, float ye, float xyint, int id)
{
	int i,k,l,mainloopexit=0,secloopexit=0;
	float xr,yr,x,y;
	int mode = 1,posn = 0, primitive=0, max_primitive=2, currentmode=0;
	float radius;
	float flux, max;
	char c, d;
	float tr[6], cont[9];
	float tempx=0., tempy=0., tempbx, tempby;

	const char *primchoice[]= {"Type : Gauss",  "Type : Rectangle", "Type : Dirac"};
	const char *modename[]={"Mode : add components", "Mode : view model"};
	tr[0]=-1.0*xyint/2.0;
	tr[1]=xyint;
	tr[2]=0.0;
	tr[3]=-1.0*xyint/2.0;
	tr[4]=0.0;
	tr[5]=xyint;

	for(i=0; i<n*n; i++)model[i] = 0.0;

	flux = 0.0;
	cpgsci(4);
	cpgslct(id);
	cpgstbg(0);
	cpgsfs(2);
	while(mainloopexit != 1)
	{
		cpgtext(-.05*xyint*n,1.05*xyint*n, "L-click/A : enter current mode, M-Click/D : change mode, R-Click/X : exit");
		cpgtext(0,xyint*n, modename[currentmode]);
		cpgcurs(&xr,&yr,&d);
		if (d=='A') { // enter current mode
			if (currentmode==0) {
				cpgsci(0);
				cpgtext(-.05*xyint*n,1.05*xyint*n, "L-click/A : enter current mode, M-Click/D : change mode, R-Click/X : exit");
				cpgsci(4);
				cpgtext(-.05*xyint*n,1.05*xyint*n, "L-click/A : add component, M-Click/D : change component type, R-Click/X : exit");
				// MODE ADD MODEL COMPONENTS
				cpgtext(0.8*xyint*n,xyint*n,primchoice[primitive]);

				while(secloopexit != 1)
				{
					cpgcurs(&xr,&yr,&c);
					/* printf("%f %f\n", xr,yr); */
					if (c == 'A') {
						/* dirac */
						if (primitive == 2) {
							cpgcirc(xr, yr, xyint);
							tempx = xr/xyint;
							tempy = yr/xyint;
							if(tempx<0.0)tempx = 0.0; if(tempy<0.0)tempy = 0.0;
							if(tempx>((float)(n-1)))tempx = (n-1);if(tempy>((float)(n-1)))tempy = (n-1);
							model[(int)tempy*n+n-1-(int)tempx]+=1.;
						} else if (primitive ==0) {
							/* gauss */
							cpgband(mode,posn,xr,yr,&x,&y,&c);
							radius = sqrt((x-xr)*(x-xr)+(y-yr)*(y-yr));
							cpgcirc(xr, yr, radius);
							cpgcirc(xr, yr, radius/2.);
							cpgcirc(xr, yr, radius/4.);
							tempx = xr/xyint;
							tempy = yr/xyint;
							if(tempx<0.0)tempx = 0.0; if(tempy<0.0)tempy = 0.0;
							if(tempx>((float)(n-1)))tempx = (n-1); if(tempy>((float)(n-1)))tempy = (n-1);
							for(k=0; k<n; k++)
							{
								for(l=0; l<n; l++)
								{
									model[k*n+n-1-l]+=exp(-((tempx-((float)l))*(tempx-((float)l)) +
									(tempy-((float)k))*(tempy-((float)k)))
									/( (radius*radius) / (xyint*xyint) ) );

								}
							}
							printf("Gaussian radius = %f\n", radius);
						} else if (primitive ==1) {
							//rectangular zone
							cpgband(mode,posn,xr,yr,&x,&y,&c);
							if (x>=xr) {tempx = x;tempbx= xr; } else { tempx = xr; tempbx= x;}
							if (y>=yr) {tempy = y;tempby= yr; } else { tempy = yr; tempby= y;}

							cpgrect(tempbx, tempx,tempby,tempy);
							tempx =tempx/xyint; tempbx =tempbx/xyint;
							tempy =tempy/xyint; tempby =tempby/xyint;
							for(k=0; k<n; k++)
							{
								for(l=0; l<n; l++)
								{
									if (((float)k>= tempby)&&((float)k< tempy)&&
										((float)l>= tempbx)&&((float)l< tempx)) model[k*n+n-1-l]+=0.5;
								}
							}

						}

					}

					if (c == 'D') {
						cpgsci(0);
						cpgtext(0.8*xyint*n,xyint*n,primchoice[primitive]);
						primitive++;
						if (primitive>max_primitive) primitive = 0;
						cpgsci(4);
						cpgtext(0.8*xyint*n,xyint*n,primchoice[primitive]);
					}
					if (c == 'X'){
						cpgsci(0);
						cpgtext(0.8*xyint*n,xyint*n,primchoice[primitive]);
						cpgsci(4);
						secloopexit = 1;
					}

				}

				secloopexit =0;


				//back to main loop, update help
				cpgsci(0);
				cpgtext(-.05*xyint*n,1.05*xyint*n,
				"L-click/A : add component, M-Click/D : change component type, R-Click/X : exit");
				cpgsci(4);
				cpgtext(-.05*xyint*n,1.05*xyint*n,
				"L-click/A : enter current mode, M-Click/D : change mode, R-Click/X : exit");

			} else if (currentmode ==1) {

				// redisplay(currentimage,n,n,xs,xe,ys,ye, 1.0, 1);
				max=0.;
				for(k=0; k<n*n; k++)
				{
					if(max<model[k])max = model[k];
				}
				cont[0] = max/100.0;
				for(k=1; k<7; k++)
				{
					cont[k] = cont[k-1]*2.;
				}
				cpgsci(4);
				cpgcont(model, n, n, 1, n, 1, n, cont, 6, tr);
			}
		}

		if (d=='D') {
			// CHANGE MODE
			cpgsci(0);
			cpgtext(0,xyint*n, modename[currentmode]);
			currentmode++;
			if (currentmode>1) currentmode=0;
			cpgsci(4);
			cpgtext(0,xyint*n, modename[currentmode]);
		}

		if (d=='X') {
			// validate model and finish
			cpgtext(0,xyint*n,"Current model validated");
			//       redisplay(currentimage,n,n,xs,xe,ys,ye,id);

			//renormalize model, display contour plot

			flux=0.;max=0.;
			for(k=0; k<n*n; k++) {
				if(model[k] <= 1e-8)model[k] = 1e-8;
				flux+=model[k];
			}
			for(k=0; k<n*n; k++)
			{
				model[k] = model[k]/flux;
				if(max<model[k])max = model[k];
			}
			cont[0] = max/100.0;
			for(k=1; k<7; k++)
			{
				cont[k] = cont[k-1]*2.;
			}
			cpgsci(4);
			cpgcont(model, n, n, 1, n, 1, n, cont, 6, tr);

			mainloopexit = 1;
		}


	}

	return SUCCESS;
}


int open_redisplay(USER *user, char* dev)
{
	user->idimage = cpgopen(dev);
	cpgslct( user->idimage );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);

	user->idv2 = cpgopen(dev);
	cpgslct( user->idv2 );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);
	
	user->idt3amp = cpgopen(dev);
	cpgslct( user->idt3amp );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);
	
	user->idt3phi = cpgopen(dev);
	cpgslct( user->idt3phi );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);

	return SUCCESS;
}


int close_redisplay(USER *user )
{
	cpgslct(user->idimage);
	cpgclos();
	cpgslct(user->idv2);
	cpgclos();
	cpgslct(user->idt3amp);
	cpgclos();
	cpgslct(user->idt3phi);
	cpgclos();
	return SUCCESS;
}

void getminmax( float* image, int nx, float *min, float *max )
{
	*min = 1e8;
	*max = -1e8;
	int i;
	for(i=0; i< nx * nx ;i++)
	{
		if( *max < image[ i ] )
			*max = image[i ];
		if( *min > image[ i ] )
			*min = image[ i ];
	}
}


int redisplay(float* image, USER *user, float displaypower, int contour)
{
	int i,j;
	float min, max;
	float tr[6];
	float C[14];
	float xyint = user->xyint;
	tr[0]= -1.0 * xyint / 2.0;
	tr[1]= xyint;
	tr[2]= 0.0;
	tr[3]= -1.0 * xyint / 2.0;
	tr[4]=0.0;
	tr[5]= xyint;

	float *img = malloc(user->iNX * user->iNY * sizeof(float));
	// Flip + power calculation
	for(i=0; i< user->iNX ;i++)
		for(j=0;j< user->iNY ; j++)
			img[i * user->iNY + user->iNX - j - 1 ] = powf( fabs(image[ i * user->iNY + j ]), displaypower);

	getminmax( img, user->iNX, &min, &max );

	cpgslct(user->idimage);
	cpgsci(1);
	cpggray(img, user->iNX, user->iNY, 1, user->iNX, 1, user->iNY, min, max, tr);
	cpgbox("BCITN",0.0,0,"BCITN",0.0,0);
	cpglab( "E-W (mas)", "S-N (mas)", "");
	if( contour > 0)
	{
		cpgsci(2);
		C[0] = max / 100.0;
		for(i=1; i<7; i++)
		{
			C[i] = C[i-1] * 2.0;
		}
		cpgcont(img, user->iNX, user->iNY, 1, user->iNX, 1, user->iNY, C, 7, tr);
	}
	free(img);


	float r, errlow, errhi;
	float complex vtemp, V0ab, V0bc,V0ca;
	char xlabel[30],ylabel[30];
	float t3max;
	
	cpgslct(user->idv2);
	cpgsci(1);
	// Powerspectrum
	if(nv2 > 0)
	{
		max = 0.0;
		min = 1e30;
		for(i=0;i<nv2;i++)
		{
			r = sqrt( u[i] * u[i] + v[i] * v[i]);
			if(max < r)max = r;
			if(min > r)min = r;
		}
		cpgenv(min, max, 0., 1.0, 0, 1);
		sprintf(xlabel,"Spatial frequency (waves)");
		sprintf(ylabel,"Powerspectrum");
		cpglab(xlabel,ylabel, "Powerspectrum");
		for(i=0;i<nv2;i++)
		{
				r = sqrt( u[i] * u[i] + v[i] * v[i] );
				cpgsci(2);
				cpgpt1( r , data[i] , 0);
				errlow = data[i] + data_err[i] ;
				errhi = data[i] - data_err[i];
				cpgerry( 1, &r, &errlow , &errhi , 1.0);
				cpgsci(4);
				cpgpt1( r , abs2( user->current_cvis[i] ) , 4);
		}
	}
	

	// CLOSURE PHASES
	if( nt3 > 0)
	  {
	    max = 0.0;
	    min = 1e30;
	    t3max = 0.;
	    for(i=0; i < nt3; i++)
	      {
		r = sqrt( u[t3in1[i]] *  u[t3in1[i]]
			  + v[t3in1[i]] *  v[t3in1[i]]
			  + u[t3in2[i]] *  u[t3in2[i]]
			  + v[t3in2[i]] *  v[t3in2[i]]
			  + u[t3in3[i]] *  u[t3in3[i]]
			  + v[t3in3[i]] *  v[t3in3[i]] );
		
		
		if(max < r)max = r;
		if(min > r)min = r;


		    if (t3max < data[nv2 +i]) t3max = data[nv2 +i];

	      }
	    
	    cpgslct(user->idt3phi);
	    cpgsci(1);
	    cpgenv(min, max, -190., 190.0, 0, 1);
	    sprintf(xlabel,"Spatial frequency (waves)");
	    sprintf(ylabel,"Closure phases");
	    cpglab(xlabel,ylabel, "Closure phases");
	    
	    
	    cpgslct(user->idt3amp);
	    cpgsci(1);
	    cpgenv(min, max, -0.001, t3max+0.001, 0, 1);
	    sprintf(xlabel,"Spatial frequency (waves)");
	    sprintf(ylabel,"Triple amplitudes");
	    cpglab(xlabel,ylabel, "Triple amplitudes");
	    
	    for(i=0;i<nt3;i++)
	      {
		    r = sqrt( u[t3in1[i]] *  u[t3in1[i]]
			      + v[t3in1[i]] *  v[t3in1[i]]
			      + u[t3in2[i]] *  u[t3in2[i]]
			      + v[t3in2[i]] *  v[t3in2[i]]
			      + u[t3in3[i]] *  u[t3in3[i]]
			      + v[t3in3[i]] *  v[t3in3[i]] );
		    
		    
		    cpgslct(user->idt3phi);
		    cpgsci(2);
		    cpgpt1( r , data[nv2 + nt3amp +i] , 0);
		    errlow = data[nv2 + nt3amp +i] + data_err[nv2 + nt3amp +i] ;
		    errhi = data[nv2 + nt3amp +i] - data_err[nv2 + nt3amp +i];
		    cpgerry( 1, &r, &errlow , &errhi , 1.0);
		    
		    if(data_err[nv2 +i] < 1e2)
		      {
			cpgslct(user->idt3amp);
			cpgsci(2);
			cpgpt1( r , data[nv2 +i] , 0);
			errlow = data[nv2 +i] + data_err[nv2 +i] ;
			errhi = data[nv2 +i] - data_err[nv2 +i];
			cpgerry( 1, &r, &errlow , &errhi , 1.0);
			
			V0ab = user->current_cvis[t3in1[i]];
			V0bc = user->current_cvis[t3in2[i]];
			V0ca = user->current_cvis[t3in3[i]];
			vtemp = V0ab * V0bc * conj(V0ca);
			
			cpgslct(user->idt3phi);
			cpgsci(4);
			cpgpt1( r , cargf( vtemp ) * 180. /M_PI , 4);
			cpgslct(user->idt3amp);
			cpgsci(4);
			cpgpt1( r , cabsf( vtemp ) , 4);
		      }

	      }
	  }
	
	return SUCCESS;
}

int dispuv( char* dev)
{
	int i;
	float max;
	char ulabel[20],vlabel[20];
	float u_point, v_point;
	max = 0.0;

	for(i=0;i<nuv;i++)
	{
		if( max < u[ i ] ) max = u[ i ];
		if( max < v[ i ] ) max = v[ i ];
	}

	sprintf( ulabel , "U /%.1e waves" , max );
	sprintf( vlabel , "V /%.1e waves" , max );

	int term;
	term= cpgopen(dev);
	cpgslct(term);
	cpgask(0);
	/* Draw UV plane coverage */
	cpgenv(-1.1, 1.1, -1.1, 1.1, 1, 1);
	cpglab(ulabel,vlabel, "UV COVERAGE");
	cpgsci(4);
	cpgsch(1.0);
	cpgpt1(-1.0,1.0,0);
	cpgptxt(-1.0,(1.0-0.02),0.0,-0.05,"v2");
	cpgsch(1.5);
	/*Draw all v2 UVs.*/
	for(i=0; i<nv2; i++)
	{
		u_point = u[i]/(max);
		v_point = v[i]/(max);
		cpgpt1(u_point,v_point,0);
		cpgpt1(-u_point,-v_point,0);
	}
	/* Draw all BS UVs */
	cpgsci(2);
	cpgsch(1.0);
	cpgpt1(-1.0,(1.0-0.07),2);
	cpgptxt(-1.0,(1.0-0.02-0.07),0.0,-0.05,"t3");
	cpgsch(1.2);
	for(i=0; i<nt3; i++)
	{
		u_point = u[t3in1[i]]/(max);
		v_point = v[t3in1[i]]/(max);
		cpgpt1(u_point,v_point,2);
		cpgpt1(-u_point,-v_point,2);
		u_point = u[t3in2[i]]/(max);
		v_point = v[t3in2[i]]/(max);
		cpgpt1(u_point,v_point,2);
		cpgpt1(-u_point,-v_point,2);
		u_point = u[t3in3[i]]/(max);
		v_point = v[t3in3[i]]/(max);
		cpgpt1(u_point,v_point,2);
		cpgpt1(-u_point,-v_point,2);
	}

	return 1;
}

//
// BSMEM additions 
//

int get_oi_fits_data(RECONST_PARAMETERS* reconst_parameters, int* status)
{
 double wavmin, wavmax, timemin, timemax;
 import_single_epoch_oifits(reconst_parameters->datafile, 1, 1, 1, 1, 1,
			    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1000, 1,
			    &wavmin, &wavmax,&timemin, &timemax);

}

int read_fits_image(char* fname, float* img, int *n, float* xyint, char* source, char* datafile, int* status)
{
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int  nfound, anynull;
	long naxes[2], fpixel, npixels;
	float nullval;

	if (*status==0)fits_open_file(&fptr, fname, READONLY, status);

	// MODIFY so that if keys do not exist, don't crash ^^

	if (*status==0)fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, status);
	//  if (*status==0)fits_read_key_str(fptr, "DATAFILE", datafile, comment, status);
	//  if (*status==0)fits_read_key_str(fptr, "TARGET", target, comment, status);
	//  if (*status==0)fits_read_key_flt(fptr, "PIXELATION", xyint, comment, status);
	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	fpixel   = 1;
	nullval  = 0;                /* don't check for null values in the image */

	//add here a compatibility check of xyint if imported as a model


	if(naxes[0] != naxes[1])
	{
		printf("Image dimension are not square.\n");
		if(*status==0)fits_close_file(fptr, status);
		return *status;
	}
	*n = naxes[0];
	/* Allocate enough memory outside of this routine */
	if(*status==0)fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval, img, &anynull, status);

	// renormalize the model !!!


	// reset the display by dispdev correctly


	if(*status==0)fits_close_file(fptr, status);

	return *status;
}


#include "oifits.c"
