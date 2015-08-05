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
 *     01/08/15 Minor tweaks for PFI
 */

#include "bsmem.h"

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
int set_memsys_dataspace( float *st , int *kb , oi_data *data , RECONST_PARAMETERS *reconst_parameters );
int reset_memsys_dataspace( float *st , int *kb , oi_data *data , RECONST_PARAMETERS *reconst_parameters );
int write_fits_image( float* img , RECONST_PARAMETERS* reconst_parameters , int* status );
void uvchuck( oi_data* data , float *st , int *kb , RECONST_PARAMETERS *reconst_parameters );
void display_oifits_info( oi_data *data );
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
  oi_data data;
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
  get_oi_fits_selection(&reconst_parameters, &status);
  get_oi_fits_data(&reconst_parameters, &data, &status);

  /**************** Map dimensions and parameters ************************************/
  user.iNX = reconst_parameters.xydim;
  user.iNY = reconst_parameters.xydim;
  user.iNUV = data.nuv;
  user.iNpow = data.npow;
  user.iNbis = data.nbis;
  user.xyint = reconst_parameters.xyint;

  /**************** Allocate MEMSYS storage areas *************************************/
  ndata = data.npow + 2 * data.nbis;
  ncells = user.iNX * user.iNY;
  nstore = 6 * ncells + 9 * ndata;
  UAREA(&ncells, &ncells, &ndata, &nstore, &level, kb);
  st = malloc(nstore * sizeof(float));

  /*************** Assign memory to work spaces *********************************/
  /* Opus/Tropus */
  user.current_visi = malloc( data.nuv * sizeof( float complex ) );
  user.current_diffvisi = malloc( data.nuv * sizeof( float complex ) );
  user.data_phasor = malloc( data.nbis * sizeof( float complex ) );

  user.filter_powerspectrum = malloc(data.npow * sizeof(int));
  user.filter_bispectrum = malloc(data.nbis * sizeof(int));
  user.UV_point = malloc(2 * data.nuv * sizeof(float));
  user.filter_UV = malloc(2 * data.nuv * sizeof(float));

  for (i = 0; i < data.npow; i++)
  {
    user.filter_powerspectrum[ i ] = 0;
  }
  for (i = 0; i < data.nuv; i++)
  {
    user.current_visi[ i ] = 0.0;
    user.current_diffvisi[ i ] = 0.0;
  }
  for (i = 0; i < data.nbis; i++)
  {
    user.filter_bispectrum[ i ] = 0;
  }
  for (i = 0; i < (2 * data .nuv); i++)
  {
    user.UV_point[ i ] = 0.0;
    user.filter_UV[ i ] = 0.0;
  }

  /*************** Initialise user.UV_point and user.bsref for transforms *********/
  user.iNUV = data.nuv;
printf("Number of unique uv points in NFFT: %d\n", user.iNUV);
  for (i = 0; i < user.iNUV; i++)
  {
    user.UV_point[ 2 * i ] = data.uv[ i ].u;
    user.UV_point[ 2 * i + 1 ] = data.uv[ i ].v;
  }
  /* make user.bsref point to data.bsref */
  user.bsref = data.bsref;

  /* Convert exp( -I*closure phase ) to Real/Imaginary parts */
  for (ii = 0; ii < data.nbis; ii++)
    user.data_phasor[ ii ] = cexp(-I * data.bisphs[ ii ] * pi / 180.0);

  /* Set up MEMSYS data area*/
  set_memsys_dataspace(st, kb, &data, &reconst_parameters);
  uvchuck(&data, st, kb, &reconst_parameters);

  int NN[ 2 ], nn[ 2 ];
  NN[ 0 ] = user.iNX;
  nn[ 0 ] = 2 * NN[ 0 ];
  NN[ 1 ] = user.iNX;
  nn[ 1 ] = 2 * NN[ 1 ];
  nfft_init_guru(&user.p, 2, NN, user.iNUV, nn, 6, PRE_FULL_PSI | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
      FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  for (i = 0; i < user.iNUV; i++)
  {
    user.p.x[ 2 * i ] = data.uv[ i ].v * user.xyint * MAS;
    user.p.x[ 2 * i + 1 ] = data.uv[ i ].u * user.xyint * MAS;
  }

  nfft_precompute_one_psi(&user.p);
  /*
  DFT_table = malloc( user.iNUV *  user.iNX * user.iNY * sizeof( float complex ) ) ;
  if( DFT_table == NULL )
    {
      printf("DFT allocation failed \n");
      getchar();
    }

  float cvfwhm = 0., bws;
  for(uu=0 ; uu < user.iNUV; uu++)
    {
      for(ii=0; ii < user.iNX; ii++)
        {
	  for(jj=0; jj < user.iNY; jj++)
	    {
	      bws = sinc( pi * MAS * user.xyint * data.uv[uu].bandwidth / data.uv[uu].wavelength
			  * ( data.uv[ uu ].u * (double)ii + data.uv[ uu ].v * (double)jj ));
	      DFT_table[ user.iNX * user.iNY * uu + ii + jj * user.iNX ] =
	        bws
		*	exp(-pi * pi * MAS * MAS * user.xyint * user.xyint / 4.0 / log(2.) * (data.uv[ uu ].u * data.uv[ uu ].u
	      			+ data.uv[ uu ].v * data.uv[ uu ].v ) * (cvfwhm * cvfwhm))
	      * cexp( - 2.0 * I * pi * MAS * user.xyint * data.uv[ uu ].u  * (double)ii )
	      * cexp( - 2.0 * I * pi * MAS * user.xyint * data.uv[ uu ].v  * (double)jj )  ;
	    }
	}
    }
  */
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
  //op_trop_check( user.iNX * user.iNY, ndata );

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

  float consistency;
  for (i = 1; i < 5; i++)
  {
    MEMTRQ(&st[ 0 ], &i, &consistency);
    if (consistency > 1e-5)
    {
      printf("Consistency error %f \n", consistency);
      //goto EXIT;
    }
  }

  /*
   // Circular object from powerspectrum
   float* circobj = malloc(  user.iNX * user.iNY * sizeof( float ) );
   // clear the coefficients
   for(i = 0 ; i < user.iNUV; i++)
   user.p.f[i] = 0.0;
   //fill those that can be with amplitudes
   for(i = 0 ; i < user.iNpow; i++)
   if( data.pow[i] > 0. )
   user.p.f[i] = sqrt( data.pow[i] );
   //do the inverse transform
   nfft_adjoint(&user.p);

   for(i=0 ; i <  user.iNX * user.iNX ; i++)
   if ( creal( user.p.f_hat[i] ) > 0.0 )
   circobj[i] = creal( user.p.f_hat[i] );
   else
   circobj[i] = 0.0;
   display( circobj, user.iNX);
   getchar();
   */
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
        dispuv(data.uv, data.bsref, data.nuv, data.npow, data.nbis, dispdev);
        break;
      case 5:
        display_command_help();
        break;
      case 6:
        if (disp_id == 0)
        {
          open_redisplay(&user, dispdev);
          redisplay(&st[ 0 ], &user, &data, reconst_parameters.displaypower, contour);
          get_default_model(startmod, &st[ 0 ], user.iNX, 0.0, user.iNX * user.xyint, 0.0, user.iNY * user.xyint, user.xyint, disp_id);
          redisplay(startmod, &user, &data, reconst_parameters.displaypower, contour);
          close_redisplay(&user);
          disp_id = 0;
        }
        else
        {
          open_redisplay(&user, dispdev);
          redisplay(&st[ 0 ], &user, &data, reconst_parameters.displaypower, contour);
          get_default_model(startmod, &st[ 0 ], user.iNX, 0.0, user.iNX * user.xyint, 0.0, user.iNY * user.xyint, user.xyint, disp_id);
          redisplay(startmod, &user, &data, reconst_parameters.displaypower, contour);
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
        redisplay(&errmap[ 0 ], &user, &data, reconst_parameters.displaypower, contour);
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
        redisplay(&st[ 0 ], &user, &data, reconst_parameters.displaypower, contour);
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
        redisplay(startmod, &user, &data, reconst_parameters.displaypower, contour);
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
          redisplay(&st[ 0 ], &user, &data, reconst_parameters.displaypower, contour);
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
            redisplay(&st[ 0 ], &user, &data, reconst_parameters.displaypower, contour);

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

        /*
         // Missing phase exploration

         int k, ntels;
         ntels = reconst_parameters.ntelescopes;
         printf( "This array uses %d telescopes \n" , ntels );
         float *missing_mat = malloc( (ntels - 1) * (ntels * (ntels - 1)) / 2
         * sizeof(float) );
         // clear matrix
         for (i = 0; i < ntels - 1; i++)
         for (j = 0; j < (ntels * (ntels - 1)) / 2; j++)
         missing_mat[ i + ntels * j ] = 0.;

         k = 0;
         for (i = 0; i < ntels - 1; i++)
         {
         for (j = i + 1; j < ntels; j++)
         {
         if (i > 0)
         missing_mat[ i - 1 + k * ntels ] = -1.;
         if (j > 0)
         missing_mat[ j - 1 + k * ntels ] = 1.;
         k++;
         }
         }

         for (j = 0; j < (ntels * (ntels - 1)) / 2; j++)
         {
         for (i = 0; i < ntels - 1; i++)
         printf( "%d \t" , (int) missing_mat[ i + ntels * j ] );
         printf( "\n" );
         }
         free( missing_mat );

         // Generate random missing phase parameters
         unsigned Rand[ 113 ];
         RanInit( Rand , 1 );
         float *missing_params = malloc( nsnapshots * (ntels - 1) * sizeof(float) );
         float *missing_phases = malloc( nsnapshots * (ntels * (ntels - 1)) / 2
         * sizeof(float) );

         // Set up gridding routine
         nfft_plan missing_plan;
         nfft_init_guru( missing_plan , 2 , NN , nsnapshots * (ntels * (ntels - 1)) / 2 ,
         nn , 6 , PRE_FULL_PSI | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT
         | FFT_OUT_OF_PLACE , FFTW_ESTIMATE | FFTW_DESTROY_INPUT );

         for (tt = 0; tt < nsnapshots; tt++)
         for (i = 0; i < (ntels * (ntels - 1)) / 2; i++)
         u = 0.;
         v = 0.;
         missing_plan.x[ 2 * i ] = v * user.xyint * MAS;
         missing_plan.x[ 2 * i + 1 ] = u * user.xyint * MAS;

         nfft_precompute_one_psi( missing_plan );

         // Get current visibilities
         // As this may not correspond to any data point, we just get them from the current image
         for (ii = 0; ii < user.iNX * user.iNX; ii++)
         missing_plan.f_hat[ ii ] = st[ ii ] + I * 0.0;
         nfft_trafo( missing_plan );
         float complex *start_visi = malloc( nsnapshots * ( ntels * (ntels - 1) ) / 2 * sizeof( float complex ) );

         // Generate missing phase parameters
         for (tt = 0; tt < nsnapshots; tt++)
         for (j = 0; j < ntels - 1; j++)
         missing_params[ tt + j * nsnapshots ] = RanGauss( Rand ) * 2. * PI;

         // Compute missing phases
         for (i = 0; i < nsnapshots * (ntels * (ntels - 1)) / 2; i++)
         missing_phases[ i ] = 0.;

         for (tt = 0; tt < nsnapshots; tt++)
         for (j = 0; j < (ntels * (ntels - 1)) / 2; j++)
         for (i = 0; i < ntels - 1; i++)
         missing_phases[ j + tt * nsnapshots ] += missing_mat[ i + ntels * j ]
         * missing_params[ i + ntels * j ];

         float complex *current_visi = malloc( data.nuv * sizeof( float complex ) );
         int uu;
         for (uu = 0; uu < nuv; uu++)
         missing_plan.f[ uu ] = start_visi[ uu ] * (cexp( I * missing_phases[ uu ] )
         - 1.);

         // Compute differential image and add it to the current image
         nfft_adjoint( missing_plan );
         for (i = 0; i < user.iNX * user.iNX; i++)
         st[ i ] += creal( missing_plan.f_hat[ i ] );

         free( missing_mat );
         free( missing_params );
         free( missing_phases );
         free( new_visi );
         getchar( );
         */
        break;

      default:
        printf("Unrecognized command\n");
        break;

    }
  } while (exit == 0);

  nfft_finalize(&user.p);
  free_oi_data(&data);
  free(st);
  free(startmod);
  free(user.current_visi);
  free(user.current_diffvisi);
  //free(user.bsref);
  free(user.filter_powerspectrum);
  free(user.filter_bispectrum);
  free(user.UV_point);
  free(user.filter_UV);
  EXIT: printf("\n\n");

  return SUCCESS;
}

int read_model( char* filename , float* mod , int N , float normalisation )
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

void uvchuck( oi_data *data , float *st , int *kb , RECONST_PARAMETERS *reconst_parameters )
{
  /* Scale UV coords, chucking amplitudes and power/bispectrum points
   * which are outside the spatial frequency cut-off set by the
   * minimum oversampling. Also discard points not wanted by the user.
   *     1990    DFB Written in Fortran.
   *     2003    HT  Translated to C and made compatible
   *             with the OI-FITS format.
   *     2005    FB  Fixed code - well, only partially...
   */

  int i;
  int ichuck, pchuck, bchuck;
  float uvmin, uvmax, uscale, vscale, ucut, vcut, ovsamp;
  float ftmp;

  /* Set initial UV selection to all UV points with all pow. and bis. */
  for (i = 0; i < user.iNUV; i++)
  {
    user.filter_UV[ 2 * i ] = 1.0;
    user.filter_UV[ 2 * i + 1 ] = 1.0;
  }
  for (i = 0; i < user.iNpow; i++)
  {
    user.filter_powerspectrum[ i ] = 1;
  }
  for (i = 0; i < user.iNbis; i++)
  {
    user.filter_bispectrum[ i ] = 1;
  }

  /* Find maximum UV coord */
  uvmax = 0.0;
  uvmin = 1e20;
  for (i = 1; i < user.iNUV; i++)
  {
    if (user.filter_UV[ i ] > 0.0)
    {
      ftmp = square(user.UV_point[ 2 * i ]) + square(user.UV_point[ 2 * i + 1 ]);
      if (ftmp > uvmax)
        uvmax = ftmp;
      if (ftmp < uvmin)
        uvmin = ftmp;
    }
  }
  uvmax = sqrt(uvmax);
  uvmin = sqrt(uvmin);
  printf("UV range:\t\t%8.0f - %8.0f wavelengths \n", uvmin, uvmax);

  /* Automatic choice of pixellation if xyint <= 0. UV cut off such that
   * the oversampling >= defined value ( OVERSAMPLING in constants.h ).
   * xyint in milliarcsec. uv in wavelengths.
   */
  printf("Array resolution:\t %f mas\n", 1. / (2.0 * MAS * uvmax));

  if (user.xyint <= 0.0)
  {
    user.xyint = 1.0 / (2.0 * MAS * OVERSAMPLING * uvmax);
    printf("Pixel size:\t\tAutomatic, %f mas\n", user.xyint);
    reconst_parameters->xyint = user.xyint;
    printf("Recommended size:\t %d pixels\n", (int) (powf(2.0, ceil(log(uvmax / uvmin * 2. * OVERSAMPLING + 1.) / log(2.) ) )));
  }
  else
    printf("Pixel size:\t\tUser defined, %f mas\n", user.xyint);

  uscale = 2.0 * (float) (user.iNX) * user.xyint * MAS;
  vscale = 2.0 * (float) (user.iNY) * user.xyint * MAS;

  ucut = (float) user.iNX / OVERSAMPLING;
  vcut = (float) user.iNY / OVERSAMPLING;

  ovsamp = 1.0 / (uvmax * user.xyint * MAS);

  printf("Image width:\t\t%d pixels, %f mas\n", user.iNX, (float) user.iNX * user.xyint);
  printf("Pix/fastest fringe:\t%f\n", ovsamp);
  /* Scale UV coords  */
  ichuck = 0;
  for (i = 0; i < user.iNUV; i++)
  {
    if (user.filter_UV[ 2 * i ] > 0.0)
    {
      user.UV_point[ 2 * i ] = uscale * user.UV_point[ 2 * i ];
      user.UV_point[ 2 * i + 1 ] = vscale * user.UV_point[ 2 * i + 1 ];

      if ((fabs(user.UV_point[ 2 * i ]) > ucut) || (fabs(user.UV_point[ 2 * i + 1 ]) > vcut))
      {
        user.filter_UV[ 2 * i ] = 0.0;
        user.filter_UV[ 2 * i + 1 ] = 0.0;
        ichuck++;
      }
    }
  }

  if (ichuck != 0)
  {
    printf("UV points tapered:\t%d ", ichuck);

    /************ MARK pow and bis if contain marked UV *********/
    pchuck = 0;
    bchuck = 0;
    /* Mark powerspectrum */
    for (i = 0; i < data->npow; i++)
    {
      if (user.filter_UV[ i ] <= 0.0)
      {
        user.filter_powerspectrum[ i ] = 0;
        pchuck++;
      }
    }
    /* Mark bispectrum */
    for (i = 0; i < data->nbis; i++)
    {
      if ((user.filter_UV[ 2 * (data ->bsref[ i ].ab.uvpnt) ] <= 0.0) || (user.filter_UV[ 2 * (data->bsref[ i ].bc.uvpnt) ] <= 0.0)
          || (user.filter_UV[ 2 * (data->bsref[ i ].ca.uvpnt) ] <= 0.0))
      {
        user.filter_bispectrum[ i ] = 0;
        bchuck++;
      }
    }
    printf("Pow: %d   Bis: %d\n", pchuck, bchuck);
  }

  // Set MEMSYS accuracies to zero for flagged or incorrect data

  for (i = 0; i < user.iNpow; i++)
  {
    if (user.filter_powerspectrum[ i ] == 0)
    {
      st[ kb[ 21 ] + i ] = 0.0;
    }
  }
  for (i = 0; i < user.iNbis; i++)
  {
    if (user.filter_bispectrum[ i ] == 0)
    {
      st[ kb[ 21 ] + user.iNpow + 2 * i ] = 0.0;
      st[ kb[ 21 ] + user.iNpow + 2 * i + 1 ] = 0.0;

      //  printf("Bispectrum %d flagged\n", i);
    }
  }

  // If requested turn off Powersp and/or Bisp points by setting accuracies to 0

  /*
   for( i=0; i<data->npow; i++ )
   {
   if( reconst_parameters->novis2 == 1 )
   {
   user.filter_powerspectrum[ i ] = 0;
   st[ kb[ 21 ] + i ]=0.0;
   }
   }

   for( i=0; i<data->nbis; i++ )
   {
   if( reconst_parameters->not3amp == 1 ) st[ kb[ 21 ] + data->npow + 2 * i ] = 0.0; // BS amp turned off
   if( reconst_parameters->not3phi == 1 ) st[ kb[ 21 ] + data->npow + 2 * i + 1 ] = 0.0; // Closures turned off
   }
   */

}

int set_memsys_dataspace( float *st , int *kb , oi_data *data , RECONST_PARAMETERS *reconst_parameters )
{
  // Note: errors bar in the data should respect the strict OIFITS definition
  // ( cf PASP paper ) - in particular, closures may be >180
  int i;
  int warn_extrapolation = 0;
  float err_pow, err_rad, err_tan, err_abs, err_phi=0;
  float pow1, powerr1, pow2, powerr2, pow3, powerr3, sqamp1, sqamp2, sqamp3, sqamperr1, sqamperr2, sqamperr3;
  printf("Loading data into memory\n");

  for (i = 0; i < data->npow; i++)
  {
    st[ kb[ 20 ] + i ] = data->pow[ i ];

    err_pow = reconst_parameters->v2a * data->powerr[ i ] + reconst_parameters->v2b;
    st[ kb[ 21 ] + i ] = 1.0 / err_pow;

    if( (data->powerr[ i ] <= 0.0 )|| isnan(data->powerr[ i ]) )
    {
      printf("Warning, error on powerspectrum %d <= 0 or NaN -- error set to infinity\n", i);
      st[ kb[ 20 ] + i ] = 0.0 ;
      st[ kb[ 21 ] + i ] = 0.0 ;
    }
  }

  for (i = 0; i < data->nbis; i++)
  {

    if  ( (reconst_parameters->forced_extrapolation == 1) || ((data->bisamperr[ i ] <= 0.) || (data->bisamperr[ i ] > (infinity - 1)))) // Missing triple amplitudes
    {
      if ((data->bsref[ i ].ab.uvpnt < data->npow) && (data->bsref[ i ].bc.uvpnt < data->npow) && (data->bsref[ i ].ca.uvpnt < data->npow))

      // if corresponding powerspectrum points are available
      {
	// Derive pseudo-triple amplitudes from powerspectrum data
        // First select the relevant powerspectra
        pow1 = data->pow[ data->bsref[ i ].ab.uvpnt ];
        powerr1 = data->powerr[ data->bsref[ i ].ab.uvpnt ];
        pow2 = data->pow[ data->bsref[ i ].bc.uvpnt ];
        powerr2 = data->powerr[ data->bsref[ i ].bc.uvpnt ];
        pow3 = data->pow[ data->bsref[ i ].ca.uvpnt ];
        powerr3 = data->powerr[ data->bsref[ i ].ca.uvpnt ];
	// printf("AB %d BC %d CA %d \n", data->bsref[ i ].ab.uvpnt, data->bsref[ i ].bc.uvpnt, data->bsref[ i ].ca.uvpnt);
        // printf("Pow1 %f Pow2 %f Pow3 %f \n", pow1, pow2, pow3);
        // printf("Powerr1 %f Powerr2 %f Powerr3 %f \n", powerr1, powerr2, powerr3);

        // Derive unbiased visibility amplitudes + noise variance
        sqamp1 = (pow1 + sqrt(square(pow1) + 2.0 * square(powerr1))) / 2.;
        sqamperr1 = 1. / (1. / sqamp1 + 2. * (3. * sqamp1 - pow1) / square(powerr1));
        sqamp2 = (pow2 + sqrt(square(pow2) + 2.0 * square(powerr2))) / 2.;
        sqamperr2 = 1. / (1. / sqamp2 + 2. * (3. * sqamp2 - pow2) / square(powerr2));
        sqamp3 = (pow3 + sqrt(square(pow3) + 2.0 * square(powerr3))) / 2.;
        sqamperr3 = 1. / (1. / sqamp3 + 2. * (3. * sqamp3 - pow3) / square(powerr3));

        // And form the triple amplitude statistics
        data->bisamp[ i ] = sqrt(sqamp1 * sqamp2 * sqamp3);
        data->bisamperr[ i ] = fabs(data->bisamp[ i ]) * sqrt(sqamperr1 / sqamp1 + sqamperr2 / sqamp2 + sqamperr3 / sqamp3);
	//   printf("Triple amplitude %d extrapolated from powerspectrum data T = %f \t E_T = %f: \r", i, data->bisamp[ i ], data->bisamperr[ i ]);
	if(warn_extrapolation == 0)
	  {
	    warn_extrapolation = 1;
	    printf("WARNING: at least one triple amplitude is missing and has been extrapolated from powerspectrum values\n");
	  }
      }

      else // missing powerspectrum points -> cannot extrapolate bispectrum
      {
        printf("WARNING: triple amplitude extrapolation from powerspectrum failed because of missing powerspectrum\n");
	data->bisamp[ i ] = 1.0;
        data->bisamperr[ i ] = infinity;
      }
    }

    if(reconst_parameters->biserrtype == 1 )
      {
	// Full elliptic approximation - Initially based on Meimon 2009 appendix E 2-3
	// See also WISARD manual

	if (i == 0)
	  printf("Bispectrum noise:\tImproved elliptic approximation \n");


	err_abs = reconst_parameters->t3ampa * data->bisamperr[ i ] + reconst_parameters->t3ampb ;
	err_phi = (reconst_parameters->t3phia * data->bisphserr[ i ] + reconst_parameters->t3phib ) * pi / 180.0 ;

	err_rad =  sqrt(
			0.5 * square( err_abs )          * square(1.+ exp(-2.*square(err_phi)))
		      + 0.5 * square( data->bisamp[ i ]) * square(1.- exp(-square(err_phi)))
		       );

	err_tan =  sqrt(
			0.5 * square( err_abs )          * square(1.- exp(-2.*square(err_phi)))
		      + 0.5 * square( data->bisamp[ i ]) * (1.- exp(-2.*square(err_phi)))
			 );

      }
    else
      {
	// Approximation - 1st order
	if (i == 0)
	  printf("Bispectrum noise:\tClassic elliptic approximation \n");
	err_rad = reconst_parameters->t3ampa * data->bisamperr[ i ] + reconst_parameters->t3ampb;
	err_tan = fabs(data->bisamp[ i ] * (reconst_parameters->t3phia * data->bisphserr[ i ] + reconst_parameters->t3phib) * pi / 180.0);
      }

    //
    // Set bispectrum errors
    //
    st[ kb[ 21 ] + data->npow + 2 * i ] = 1.0 / err_rad;
    st[ kb[ 21 ] + data->npow + 1 + 2 * i ] = 1.0 / err_tan;

    // if bisamperr < 0 / bisphserr < 0 then the point had been flagged, so reset accuracy to 0
    if( ( (data->bisamperr[ i ] <= 0.0) || isnan(data->bisamperr[ i ]) ) || ( (data->bisphserr[ i ] <= 0.0 ) || isnan(data->bisphserr[ i ] ) ))
      {
	// printf("Warning, error on bispectrum %d <= 0 -- error set to infinity\n", i);
	st[ kb[ 21 ] + data->npow + 2 * i ] = 0.0 ;
	st[ kb[ 21 ] + data->npow + 1 + 2 * i ] = 0.0 ;
      }

    //
    // Set Bispectrum data, rotated by exp( -i*closure_data )
    //
    if(reconst_parameters->biserrtype == 1 )
	// Improved elliptic approximation
      st[ kb[ 20 ]+data->npow + 2 * i ] = data->bisamp[ i ] * (2.- exp( - 0.5*square(err_phi) ))  ;
   else
	// Approximation 1st order
        st[ kb[ 20 ] + data->npow + 2 * i ] = data->bisamp[ i ];

    st[ kb[ 20 ] + data->npow + 1 + 2 * i ] = 0.0;

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

void display_oifits_info( oi_data *data )
{
  float ratio, nmin1 = 0.0, nmin2 = 0.0, nmin3 = 0.0, nmax1 = 1e10, nmax2 = 1e10, nmax3 = 1e10, tot1 = 0., tot2 = 0., tot3 = 0.;
  char choice[ 5 ], tempstr;
  int ii;
  printf("Display powerspectrum data ? y/[ N ] ");
  if(fgets(choice, sizeof(choice), stdin) ==NULL)
    printf("Error getting answer\n");
  tempstr = choice[ 0 ];
  for (ii = 0; ii < data->npow; ii++)
  {
    ratio = fabs(data->pow[ ii ] / data->powerr[ ii ]);
    if ((tempstr == 'y') || (tempstr == 'Y'))
      printf("N: %d Amp: %9f +/-%9f | S/N %3.1f\n", ii, data->pow[ ii ], data->powerr[ ii ], ratio);
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
	tot1 += ratio / data->npow;
      }
  }
  printf("Display bispectrum data ? y/[ N ] ");
  if(fgets(choice, sizeof(choice), stdin) ==NULL)
    printf("Error getting answer\n");

  tempstr = choice[ 0 ];
  for (ii = 0; ii < data->nbis; ii++)
  {
    ratio = fabs(data->bisamp[ ii ] / data->bisamperr[ ii ]);
    if ((tempstr == 'y') || (tempstr == 'Y'))
      printf("N: %d Amp : %9f +/-%9f | S/N %3.1f | Phs : %9f +/-%9f\n", ii, data->bisamp[ ii ], data->bisamperr[ ii ], ratio,
          data->bisphs[ ii ], data->bisphserr[ ii ]);
    if (ii == 0)
    {
      nmin2 = ratio;
      nmax2 = ratio;
      nmin3 = data->bisphserr[ ii ];
      nmax3 = data->bisphserr[ ii ];
    }
    if (ratio < nmin2)
      nmin2 = ratio;
    if (ratio > nmax2)
      nmax2 = ratio;
    if (data->bisphserr[ ii ] < nmin3)
      nmin3 = data->bisphserr[ ii ];
    if (data->bisphserr[ ii ] > nmax3)
      nmax3 = data->bisphserr[ ii ];
    tot2 += ratio / data->nbis;
    tot3 += data->bisphserr[ ii ] / data->nbis;
  }

  printf("SNR |   pow    |  t3amp    | t3phs    |\n");
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

void __stdcall VMEMEX( float *image , float *data )
{
  // Visible-to-Data transform
  // Takes present model image and calculates the appropriate non-linear
  // mock data values. (Powerpectrum and bispectrum points)

  int i;
  int namp, nbs, nuv;
  float complex vtemp, V0ab, V0bc,V0ca;
  namp = user.iNpow;
  nbs = user.iNbis;
  nuv = user.iNUV;

  user.imflux = 0.0;
  int ii, uu;
    for (ii = 0; ii < user.iNX * user.iNX; ii++)
   {
    user.p.f_hat[ ii ] = image[ ii ] + I * 0.0;
    user.imflux += image[ ii ];
   }
   nfft_trafo(&user.p);

  /*
  // DFT

  for (uu = 0; uu < nuv; uu++)
    {
      user.p.f[ uu ] = 0.;
      for (ii = 0; ii < user.iNX * user.iNX; ii++)
	user.p.f[ uu ] += DFT_table[  uu * user.iNX * user.iNX + ii] *  image[ ii ]  ;
    }

  // END DFT
  */
  for (uu = 0; uu < nuv; uu++)
    user.current_visi[ uu ] = user.p.f[ uu ];

  // printf("Check - visi 0 == %f + i * %f\n", creal(user.current_visi[10]), cimag(user.current_visi[10]));

  // Powerpesctrum
  for (i = 0; i < namp; i++)
  {
    if (user.filter_powerspectrum[ i ] > 0)
      data[ i ] = abs2(user.current_visi[ i ]);
  else
      {
	data[ i ] = 0.;
	//	printf( "WARNING - Flagged powerspectrum\n" );
      }
  }
  //printf("Check - powerspectrum 0 == %f\n", data[0]);
  // Bispectrum
  for (i = 0; i < nbs; i++)
  {

    if (user.filter_bispectrum[ i ] > 0)
    {
      V0ab = user.current_visi[ user.bsref[ i ].ab.uvpnt ];
      V0bc = user.current_visi[ user.bsref[ i ].bc.uvpnt ];
      V0ca = user.current_visi[ user.bsref[ i ].ca.uvpnt ];
      if (user.bsref[ i ].ab.sign < 0)
        V0ab = conj(V0ab);
      if (user.bsref[ i ].bc.sign < 0)
        V0bc = conj(V0bc);
      if (user.bsref[ i ].ca.sign < 0)
        V0ca = conj(V0ca);

      vtemp = V0ab * V0bc * V0ca * user.data_phasor[ i ];

      data[ namp + 2 * i ] = creal(vtemp);
      data[ namp + 2 * i + 1 ] = cimag(vtemp);

    }
    else
    {
      //  			printf( "WARNING - Flagged bispectrum\n" );
      data[ namp + 2 * i ] = 0.0;
      data[ namp + 2 * i + 1 ] = 0.0;
    }

  }

}

void __stdcall VOPUS( float *dh , float *dd )
{
  // Visible-to-Data differential transform
  // Note current_diffvisi is dVisibility, dh is dImage, dd is dData
  int i;
  int namp, nbs, nuv;
  float complex vtemp;
  float complex V0ab, V0bc, V0ca;
  float complex VISab, VISbc, VISca;
  namp = user.iNpow;
  nbs = user.iNbis;
  nuv = user.iNUV;

  int ii, uu;
  for (ii = 0; ii < user.iNX * user.iNX; ii++)
  {
    user.p.f_hat[ ii ] = dh[ ii ] + I * 0.0;
  }

  nfft_trafo(&user.p);

  /*
 for (uu = 0; uu < nuv; uu++)
    {
      user.p.f[ uu ] = 0.;
      for (ii = 0; ii < user.iNX * user.iNX; ii++)
	 user.p.f[ uu ] += DFT_table[  uu * user.iNX * user.iNX + ii] *  dh[ ii ] ;
    }
  */

  for (uu = 0; uu < nuv; uu++)
        user.current_diffvisi[ uu ] = user.p.f[ uu ];

  /* Powerspectrum */
  for (i = 0; i < namp; i++)
  {
    if (user.filter_powerspectrum[ i ] > 0)
      dd[ i ] = 2.0 * (creal(user.current_diffvisi[ i ]) * creal(user.current_visi[ i ]) + cimag(user.current_diffvisi[ i ]) * cimag(
          user.current_visi[ i ]));
  }

  /* Bispectrum */
  for (i = 0; i < nbs; i++)
  {
    if (user.filter_bispectrum[ i ] > 0)
    {
      V0ab = user.current_visi[ user.bsref[ i ].ab.uvpnt ];
      V0bc = user.current_visi[ user.bsref[ i ].bc.uvpnt ];
      V0ca = user.current_visi[ user.bsref[ i ].ca.uvpnt ];
      VISab = user.current_diffvisi[ user.bsref[ i ].ab.uvpnt ];
      VISbc = user.current_diffvisi[ user.bsref[ i ].bc.uvpnt ];
      VISca = user.current_diffvisi[ user.bsref[ i ].ca.uvpnt ];

      /* conjugate if baseline in opposite part of uv-plane. */
      if (user.bsref[ i ].ab.sign < 0)
        V0ab = conj(V0ab);
      if (user.bsref[ i ].bc.sign < 0)
        V0bc = conj(V0bc);
      if (user.bsref[ i ].ca.sign < 0)
        V0ca = conj(V0ca);
      if (user.bsref[ i ].ab.sign < 0)
        VISab = conj(VISab);
      if (user.bsref[ i ].bc.sign < 0)
        VISbc = conj(VISbc);
      if (user.bsref[ i ].ca.sign < 0)
        VISca = conj(VISca);

      /* differential response calculation */
      vtemp = user.data_phasor[ i ] * (VISab * V0bc * V0ca + VISbc * V0ca * V0ab + VISca * V0ab * V0bc);
      dd[ namp + 2 * i ] = creal(vtemp);
      dd[ namp + 2 * i + 1 ] = cimag(vtemp);

    }
  }

}

void __stdcall VTROP( float *dh , float *dd )
{ // Data-to-Visible differential transform
  int i;
  int namp, nbs, nuv;
  int ab_pnt, bc_pnt, ca_pnt;
  float complex V0ab,V0bc,V0ca;
  float complex VISab,VISbc,VISca;
  float complex bs;
  namp = user.iNpow;
  nbs = user.iNbis;
  nuv = user.iNUV;

  for (i = 0; i < nuv; i++)
    user.current_diffvisi[ i ] = 0.0;

  /* Powerspectrum contribution */
  for (i = 0; i < namp; i++)
  {
    if (user.filter_powerspectrum[ i ] > 0)
      user.current_diffvisi[ i ] += user.current_visi[ i ] * 2.0 * dd[ i ];
  }

  /* Bispectrum contribution */
  for (i = 0; i < nbs; i++)
  {
    if (user.filter_bispectrum[ i ] > 0)
    {

      /* set uv reference for simplification */
      ab_pnt = user.bsref[ i ].ab.uvpnt;
      bc_pnt = user.bsref[ i ].bc.uvpnt;
      ca_pnt = user.bsref[ i ].ca.uvpnt;

      /* previous visibilities */
      V0ab = conj(user.current_visi[ ab_pnt ]);
      V0bc = conj(user.current_visi[ bc_pnt ]);
      V0ca = conj(user.current_visi[ ca_pnt ]);

      /* conjugate if baseline in opposite part of uv-plane. */
      if (user.bsref[ i ].ab.sign < 0)
        V0ab = conj(V0ab);
      if (user.bsref[ i ].bc.sign < 0)
        V0bc = conj(V0bc);
      if (user.bsref[ i ].ca.sign < 0)
        V0ca = conj(V0ca);

      /* input bispectrum differential */
      bs = (dd[ namp + 2 * i ] + I * dd[ namp + 2 * i + 1 ]) * conj(user.data_phasor[ i ]);

      /* differential response calculation */
      VISab = bs * V0bc * V0ca;
      VISbc = bs * V0ca * V0ab;
      VISca = bs * V0ab * V0bc;

      /* conjugate if baseline in opposite part of uv-plane. */
      if (user.bsref[ i ].ab.sign < 0)
        VISab = conj(VISab);
      if (user.bsref[ i ].bc.sign < 0)
        VISbc = conj(VISbc);
      if (user.bsref[ i ].ca.sign < 0)
        VISca = conj(VISca);

      /* Add to differential response */
      user.current_diffvisi[ ab_pnt ] += VISab;
      user.current_diffvisi[ bc_pnt ] += VISbc;
      user.current_diffvisi[ ca_pnt ] += VISca;
    }
  }

  int ii, uu;
   for (uu = 0; uu < nuv; uu++)
    user.p.f[ uu ] = user.current_diffvisi[ uu ];

   nfft_adjoint(&user.p);

  /*
  for (ii = 0; ii < user.iNX * user.iNX; ii++)
    {
      user.p.f_hat[ ii ] = 0. ;
      for (uu = 0; uu < nuv; uu++)
	user.p.f_hat[ ii ] += 1. / ( DFT_table[  uu * user.iNX * user.iNX + ii] ) *  user.current_diffvisi[ uu ] ;
    }
  */
  for (ii = 0; ii < user.iNX * user.iNX; ii++)
  {
    dh[ ii ] = creal(user.p.f_hat[ ii ]);
  }

}

void op_trop_check( int npix , int ndata )
{
  int i, jj;
  float *hid1, *vis1, *hid2, *vis2, *temp;
  float sumh, sumv;
  float precision;
  unsigned Rand[ 113 ];
  RanInit(Rand, 1);// -1 for time based, +1 for non time based
  printf("Testing Opus/tropus consistency with %d pixels and %d data\n", npix, ndata);

  hid1 = malloc(npix * sizeof(float));
  vis1 = malloc(ndata * sizeof(float));
  hid2 = malloc(npix * sizeof(float));
  vis2 = malloc(ndata * sizeof(float));

  temp = malloc(ndata * sizeof(float));

  for (jj = 0; jj < 10; jj++)
  {

    for (i = 0; i < npix; i++)
    {
      hid1[ i ] = Ranfloat(Rand);
      hid2[ i ] = 0.0;
    }
    for (i = 0; i < ndata; i++)
    {
      vis1[ i ] = Ranfloat(Rand);
      vis2[ i ] = 0.0;
    }

    VMEMEX(hid1, temp);
    VOPUS(hid1, vis2);
    VTROP(hid2, vis1);
    sumh = 0.0;
    sumv = 0.0;
    for (i = 0; i < npix; i++)
      sumh += hid1[ i ] * hid2[ i ];
    for (i = 0; i < ndata; i++)
      sumv += vis1[ i ] * vis2[ i ];
    precision = fabs((sumh - sumv) / sumh);
    printf("OPTROP Test %d Precision : %e ( SUMH : %f SUMV : %f ) \n", jj, precision, sumh, sumv);
  }
 free(vis1);
 free(vis2);
 free(hid1);
 free(hid2);
 free(temp);


}

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
