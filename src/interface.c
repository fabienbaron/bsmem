/*
* interface.c
*
* The Contents of this file are made available subject to the terms
* of the GNU Lesser General Public License:
*
* Copyright (C) 2008-2009 Fabien Baron
*
* This library is free software: you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library.  If not, see
* http://www.gnu.org/licenses/
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <ctype.h>
#include "bsmem.h"
#include "cpgplot.h"

int s_equals(char *A, char *B);
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
  printf("\n\n**********  BSMEM v1.5   ******************\n");
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
      else if(strcmp(argv[i],"-v2a") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->v2a);
      else if(strcmp(argv[i],"-v2b") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->v2b);
      else if(strcmp(argv[i],"-t3ampa") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->t3ampa);
      else if(strcmp(argv[i],"-t3ampb") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->t3ampb);
      else if(strcmp(argv[i],"-t3phia") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->t3phia);
      else if(strcmp(argv[i],"-t3phib") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->t3phib);
      else if(strcmp(argv[i],"-t3phib") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->t3phib);
      else if (strcmp(argv[i],"-berr") == 0)
	sscanf(argv[i+1],"%d", &reconst_parameters->biserrtype);
      else if (strcmp(argv[i],"-ferr") == 0)
	sscanf(argv[i+1], "%f", &reconst_parameters->fluxerr);
      else if(strcmp(argv[i],"-novis2") == 0)
	reconst_parameters->novis2 = 1;
      else if(strcmp(argv[i],"-not3amp") == 0)
	reconst_parameters->not3amp = 1;
      else if(strcmp(argv[i],"-noui") == 0)
	reconst_parameters->noui = 1;
      else if(strcmp(argv[i],"-forcext") == 0)
	reconst_parameters->forced_extrapolation = 1;
      else if(strcmp(argv[i],"-not3phi") == 0)
	    reconst_parameters->not3phi = 1;
      else if(strcmp(argv[i],"-ver") == 0)
	{
	  printf("***********************************************\n");
	  printf("*          B S M E M  v1.5                    *\n");
	  printf("*                                             *\n");
	  printf("* Current version  : Fabien Baron             *\n");
	  printf("* OIFITS library   : John Young               *\n");
	  printf("* Based on code by : David Buscher            *\n");
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
      printf("-v2a:\t\t Multiplicative scaling factor for powerspectrum errors (e'= a * e + b).\n");
      printf("-v2b:\t\t Additive factor for powerspectrum errors (e'= a * e + b).\n");
      printf("-t3ampa:\t Multiplicative scaling factor for triple amplitudes (e'= a * e + b).\n");
      printf("-t3ampb:\t Additive factor for triple amplitudes (e'= a * e + b).\n");
      printf("-t3phia:\t Multiplicative scaling factor for closure phases (e'= a * e + b).\n");
      printf("-t3phib:\t Additive factor for closure phases (e'= a * e + b).\n");
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
	reconst_parameters->bisio = 0;
	reconst_parameters->powio = 0;
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
	reconst_parameters->v2a = 1.0;
	reconst_parameters->t3ampa = 1.0;
	reconst_parameters->t3phia = 1.0;
	reconst_parameters->v2b = 0.0;
	reconst_parameters->t3ampb = 0.0;
	reconst_parameters->t3phib = 0.0;
	reconst_parameters->biserrtype = 0;
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

	user->idpow = cpgopen(dev);
	cpgslct( user->idpow );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);
	
	user->idt3amp = cpgopen(dev);
	cpgslct( user->idt3amp );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);
	
	user->idt3phs = cpgopen(dev);
	cpgslct( user->idt3phs );
	cpgask(0);
	cpgenv( 0.0 , user->iNX * user->xyint, 0.0 , user->iNY * user->xyint, 1, -2);

	return SUCCESS;
}


int close_redisplay(USER *user )
{
	cpgslct(user->idimage);
	cpgclos();
	cpgslct(user->idpow);
	cpgclos();
	cpgslct(user->idt3amp);
	cpgclos();
	cpgslct(user->idt3phs);
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


int redisplay(float* image, USER *user, oi_data *data, float displaypower, int contour)
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
	
	cpgslct(user->idpow);
	cpgsci(1);
	// Powerspectrum
	if(data->npow > 0)
	{
		max = 0.0;
		min = 1e30;
		for(i=0;i<data->npow;i++)
		{
			r = sqrt( data->uv[i].u * data->uv[i].u + data->uv[i].v * data->uv[i].v );
			if(max < r)max = r;
			if(min > r)min = r;
		}
		cpgenv(min, max, 0., 1.0, 0, 1);
		sprintf(xlabel,"Spatial frequency (waves)");
		sprintf(ylabel,"Powerspectrum");
		cpglab(xlabel,ylabel, "Powerspectrum");
		for(i=0;i<data->npow;i++)
		{
			if ( user->filter_powerspectrum[i] > 0 )
			{
				r = sqrt( data->uv[i].u * data->uv[i].u + data->uv[i].v * data->uv[i].v );
				cpgsci(2);
				cpgpt1( r , data->pow[i] , 0);
				errlow = data->pow[i] + data->powerr[i] ;
				errhi = data->pow[i] - data->powerr[i];
				cpgerry( 1, &r, &errlow , &errhi , 1.0);
				cpgsci(4);
				cpgpt1( r , abs2( user->current_visi[i] ) , 4);
			}
		}
	}
	

	// CLOSURE PHASES
	if( data->nbis > 0)
	  {
	    max = 0.0;
	    min = 1e30;
	    t3max = 0.;
	    for(i=0; i < data->nbis; i++)
	      {
		r = sqrt( data->uv[data->bsref[i].ab.uvpnt].u *  data->uv[data->bsref[i].ab.uvpnt].u
			  + data->uv[data->bsref[i].ab.uvpnt].v *  data->uv[data->bsref[i].ab.uvpnt].v
			  + data->uv[data->bsref[i].bc.uvpnt].u *  data->uv[data->bsref[i].bc.uvpnt].u
			  + data->uv[data->bsref[i].bc.uvpnt].v *  data->uv[data->bsref[i].bc.uvpnt].v
			  + data->uv[data->bsref[i].ca.uvpnt].u *  data->uv[data->bsref[i].ca.uvpnt].u
			  + data->uv[data->bsref[i].ca.uvpnt].v *  data->uv[data->bsref[i].ca.uvpnt].v );
		
		
		if(max < r)max = r;
		if(min > r)min = r;

		if(user->filter_bispectrum[i]>0)
		  {
		    if (t3max < data->bisamp[i]) t3max = data->bisamp[i];
		  }		

	      }
	    
	    cpgslct(user->idt3phs);
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
	    
	    for(i=0;i<data->nbis;i++)
	      {
		if(user->filter_bispectrum[i]>0)
		  {
		    r = sqrt( data->uv[data->bsref[i].ab.uvpnt].u *  data->uv[data->bsref[i].ab.uvpnt].u
			      + data->uv[data->bsref[i].ab.uvpnt].v *  data->uv[data->bsref[i].ab.uvpnt].v
			      + data->uv[data->bsref[i].bc.uvpnt].u *  data->uv[data->bsref[i].bc.uvpnt].u
			      + data->uv[data->bsref[i].bc.uvpnt].v *  data->uv[data->bsref[i].bc.uvpnt].v
			      + data->uv[data->bsref[i].ca.uvpnt].u *  data->uv[data->bsref[i].ca.uvpnt].u
			      + data->uv[data->bsref[i].ca.uvpnt].v *  data->uv[data->bsref[i].ca.uvpnt].v );
		    
		    
		    cpgslct(user->idt3phs);
		    cpgsci(2);
		    cpgpt1( r , data->bisphs[i] , 0);
		    errlow = data->bisphs[i] + data->bisphserr[i] ;
		    errhi = data->bisphs[i] - data->bisphserr[i];
		    cpgerry( 1, &r, &errlow , &errhi , 1.0);
		    
		    if(data->bisamperr[i] < 1e2)
		      {
			cpgslct(user->idt3amp);
			cpgsci(2);
			cpgpt1( r , data->bisamp[i] , 0);
			errlow = data->bisamp[i] + data->bisamperr[i] ;
			errhi = data->bisamp[i] - data->bisamperr[i];
			cpgerry( 1, &r, &errlow , &errhi , 1.0);
			
			V0ab = user->current_visi[user->bsref[i].ab.uvpnt];
			V0bc = user->current_visi[user->bsref[i].bc.uvpnt];
			V0ca = user->current_visi[user->bsref[i].ca.uvpnt];
			if(user->bsref[i].ab.sign<0)V0ab = conj(V0ab);
			if(user->bsref[i].bc.sign<0)V0bc = conj(V0bc);
			if(user->bsref[i].ca.sign<0)V0ca = conj(V0ca);
			vtemp = V0ab * V0bc * V0ca;
			
			cpgslct(user->idt3phs);
			cpgsci(4);
			cpgpt1( r , cargf( vtemp ) * 180. /pi , 4);
			cpgslct(user->idt3amp);
			cpgsci(4);
			cpgpt1( r , cabsf( vtemp ) , 4);
		      }
		  }
	      }
	  }
	
	return SUCCESS;
}

int dispuv(oi_uv* A, oi_bsref* B, int nuv, int npow, int nbis, char* dev)
{
	int i;
	float max;
	float u,v;
	char ulabel[20],vlabel[20];

	max = 0.0;

	for(i=0;i<nuv;i++)
	{
		if( max < A[ i ].u ) max = A[ i ].u;
		if( max < A[ i ].v ) max = A[ i ].v;
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
	cpgptxt(-1.0,(1.0-0.02),0.0,-0.05,"powerspectrum");
	cpgsch(1.5);
	/*Draw all pow UVs.*/
	for(i=0; i<npow; i++)
	{
		u = A[i].u/(max);
		v = A[i].v/(max);
		cpgpt1(u,v,0);
		cpgpt1(-u,-v,0);
	}
	/* Draw all BS UVs */
	cpgsci(2);
	cpgsch(1.0);
	cpgpt1(-1.0,(1.0-0.07),2);
	cpgptxt(-1.0,(1.0-0.02-0.07),0.0,-0.05,"bispectrum");
	cpgsch(1.2);
	for(i=0; i<nbis; i++)
	{
		u = A[B[i].ab.uvpnt].u/(max);
		v = A[B[i].ab.uvpnt].v/(max);
		cpgpt1(u,v,2);
		cpgpt1(-u,-v,2);
		u = A[B[i].bc.uvpnt].u/(max);
		v = A[B[i].bc.uvpnt].v/(max);
		cpgpt1(u,v,2);
		cpgpt1(-u,-v,2);
		u = A[B[i].ca.uvpnt].u/(max);
		v = A[B[i].ca.uvpnt].v/(max);
		cpgpt1(u,v,2);
		cpgpt1(-u,-v,2);
	}

	return 1;
}
