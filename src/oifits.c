/*
* oifits.c
*
* The Contents of this file are made available subject to the terms
* of the GNU Lesser General Public License:
*
* Copyright (C) 2003-2010 Fabien Baron, Hrobjartur Thorsteinsson
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

/* Oifits routines.
* This package uses the Oifits Exchange routines by John Young to view,
* select and extract oi data.
*/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bsmem.h"

#define uv_threshold 5.0e-5

int compare_uv(oi_uv uv, oi_uv withuv, float thresh)
{
	float sqthresh;
	float pcompu,pcompv;
	float mcompu,mcompv;


	sqthresh = thresh*thresh;

	pcompu = 2.0*(uv.u-withuv.u)*(uv.u-withuv.u)/(uv.u*uv.u+withuv.u*withuv.u);
	pcompv = 2.0*(uv.v-withuv.v)*(uv.v-withuv.v)/(uv.v*uv.v+withuv.v*withuv.v);
	mcompu = 2.0*(uv.u+withuv.u)*(uv.u+withuv.u)/(uv.u*uv.u+withuv.u*withuv.u);
	mcompv = 2.0*(uv.v+withuv.v)*(uv.v+withuv.v)/(uv.v*uv.v+withuv.v*withuv.v);

	/* To handle zeros */
	if((uv.u-withuv.u) == 0.0)
	{
		pcompu = 0.0;
	}
	if((uv.v-withuv.v) == 0.0)
	{
		pcompv = 0.0;
	}
	if((uv.u+withuv.u) == 0.0)
	{
		mcompu = 0.0;
	}
	if((uv.v+withuv.v) == 0.0)
	{
		mcompv = 0.0;
	}


	/* If uv is same as withuv then return 1 */
	if((pcompu<sqthresh)&&(pcompv<sqthresh))
	{
		return 1;
	}
	/* If uv is same as withuv but conjugated then return -1 */
	if((mcompu<sqthresh)&&(mcompv<sqthresh))
	{
		return -1;
	}

	return 0;
}


int get_oi_fits_data(RECONST_PARAMETERS* reconst_parameters, oi_data* data, int* status)
{
	/*
	*  ROUTINE FOR READING OI_FITS DATA.
	*  All uv coords are conjugated into +ve u half plane.
	*/
	/* Declare data structure variables*/
	/*     OI-FITS */
	/* oi_array array; */
	oi_wavelength wave;
	oi_vis2 vis2;
	oi_t3 t3;
	/* Declare other variables */
	int i,k,l;
	int phu;
	int npve = 0, uvexists;
	oi_uv puv1,puv2,puv3;
	fitsfile *fptr;

	/* If error do nothing */
	if(*status) return *status;

	/* Read fits file */
	fits_open_file(&fptr, reconst_parameters->datafile, READONLY, status);
	if(*status) {
		fits_report_error(stderr, *status);
		exit(1);
	}

	/* Allocate memory */
	data->pow = malloc( (reconst_parameters->numvis2 + 1) *sizeof(float));
	data->powerr = malloc( (reconst_parameters->numvis2 + 1 ) *sizeof(float));
	data->bisamp = malloc(reconst_parameters->numt3 *sizeof(float));
	data->bisamperr = malloc(reconst_parameters->numt3 *sizeof(float));
	data->bisphs = malloc(reconst_parameters->numt3 *sizeof(float));
	data->bisphserr = malloc(reconst_parameters->numt3 *sizeof(float));
	data->bsref = malloc(reconst_parameters->numt3 *sizeof(oi_bsref));
	data->uv = malloc((1 + reconst_parameters->numvis2+3 *reconst_parameters->numt3)*sizeof(oi_uv));
	/* Allocate as much space for UV as possible initially and then reallocate in the end */

	/* Read in visibility */
	if(*status==0)
	{
		/* powerspectrum */
		fits_movabs_hdu(fptr,1,NULL,status);

		/* Add total "zero" flux component for normalization */
		data->npow = 1;
		data->nuv  = 1;
		data->pow[0]=1.0;
		data->powerr[0]=reconst_parameters->fluxerr;
		data->uv[0].u=0.;
		data->uv[0].v=0.;
		while(*status==0)
		{
			read_next_oi_vis2(fptr, &vis2, status);
			fits_get_hdu_num(fptr, &phu);
			read_oi_wavelength(fptr, vis2.insname, &wave, status);
			fits_movabs_hdu(fptr,phu,NULL,status);

			if(*status==0)
			{
				if(vis2.record[0].target_id == reconst_parameters->target_id)
				{
					for(i=0; i<vis2.numrec; i++)
					{
					  printf("Reading V2 record %d/%ld\r", i+1, vis2.numrec);
						for(k=0; k<vis2.nwave; k++)
						{
							if(((wave.eff_wave[k]*billion)> reconst_parameters->minband)
									&&((wave.eff_wave[k]*billion)< reconst_parameters->maxband)
							   &&(!(vis2.record[i].flag[k])))
							{
							  
							      data->pow[data->npow] = (float)vis2.record[i].vis2data[k];
							      data->powerr[data->npow] = (float)vis2.record[i].vis2err[k];
							      
							      data->uv[data->nuv].u = (float)(vis2.record[i].ucoord / wave.eff_wave[k]);
							      data->uv[data->nuv].v = (float)(vis2.record[i].vcoord / wave.eff_wave[k]);
							      //printf("%d %e %e %e\n", data->nuv,  wave.eff_wave[k], wave.eff_band[k], wave.eff_band[k] / wave.eff_wave[k]);
							      data->uv[data->nuv].wavelength = wave.eff_wave[k];
							      data->uv[data->nuv].bandwidth = wave.eff_band[k];
							      
							      /* flip into +u half plane */
							      if(data->uv[data->nuv].u<0.0)
								{
								  data->uv[data->nuv].u = -data->uv[data->nuv].u;
								  data->uv[data->nuv].v = -data->uv[data->nuv].v;
								}
							  
							      data->npow++;
							      data->nuv++;
						
							}
						}
					}
					printf("\n");
				}
			}
			/* free memory */
			if(*status == 0)
			{
				free_oi_wavelength(&wave);
				free_oi_vis2(&vis2);
			}
		}
		*status = 0;

		/* bispectrum */
		fits_movabs_hdu(fptr,1,NULL,status);
		data->nbis  = 0;
		while(*status==0)
		{
			read_next_oi_t3(fptr, &t3, status);
			fits_get_hdu_num(fptr, &phu);
			read_oi_wavelength(fptr, t3.insname, &wave, status);
			fits_movabs_hdu(fptr,phu,NULL,status);

			if(*status==0)
			{
				if(t3.record[0].target_id == reconst_parameters->target_id)
				{
					for(i=0; i<t3.numrec; i++)
					{
					printf("Reading T3 record %d/%ld\r", i+1, t3.numrec);
					  for(k=0; k<t3.nwave; k++)
						{
						  //	  printf("i: %d k: %d\t ", i, k); 
							if(((wave.eff_wave[k]*billion)>reconst_parameters->minband)
							   &&((wave.eff_wave[k]*billion)<reconst_parameters->maxband)
							   &&(!(t3.record[i].flag[k])))
							  {
							    /* Trick to use closure data without available bis amplitudes */
							    if(isnan(t3.record[i].t3amp[k]))
								{
									data->bisamp[data->nbis] = 1.0;
									data->bisamperr[data->nbis] = infinity ;
								}
								else
								{
									data->bisamp[data->nbis] = (float)(t3.record[i].t3amp[k]);
									data->bisamperr[data->nbis] = (float)(t3.record[i].t3amperr[k]);
								}

							    // Trick to use amplitudes without any closures
							    if(isnan(t3.record[i].t3phi[k]))
							    {
									data->bisphs[data->nbis] = 0.0;
									data->bisphserr[data->nbis] = infinity;
							    }
							    else
							    {
									data->bisphs[data->nbis] = (float)(t3.record[i].t3phi[k]);
									data->bisphserr[data->nbis] = (float)(t3.record[i].t3phierr[k]);
							    }

								// data->bistime[data->nbis] = t3.time;

								/* Read UV coords and check if they exist. If do not exist -> update UV.
								* Set the bsref.
								*/
								puv1.u = (float)( t3.record[ i ].u1coord / wave.eff_wave[ k ] );
								puv1.v = (float)( t3.record[ i ].v1coord / wave.eff_wave[ k ] );
								puv2.u = (float)( t3.record[ i ].u2coord / wave.eff_wave[ k ] );
								puv2.v = (float)( t3.record[ i ].v2coord / wave.eff_wave[ k ] );
								puv3.u = -(puv1.u + puv2.u);
								puv3.v = -(puv1.v + puv2.v);
								//	if(i == 15) 
								//  printf("%f %f %f %f %f %f %f \n",puv1.u, puv1.v, puv2.u, puv2.v, puv3.u, puv3.v, wave.eff_wave[ k ] );
								/* Check if UV1, UV2, UV3 exist */

								/*uv1*/
								uvexists = 0;
								for(l=0; l<data->nuv; l++)
								  {
								    /*if((l==361)&&(i==15))
								      printf("\n k: %d uvexist %d l %d uv.u %f uv.v %f uvl.u %f uvl.v %f dist %f\n", k, uvexists,
								    	     l, puv1.u, puv1.v,data->uv[l].u,data->uv[l].v,
								    	     sqrt( (-puv1.u-data->uv[l].u)*(-puv1.u-data->uv[l].u) + (-puv1.v-data->uv[l].v)*(-puv1.v-data->uv[l].v) ) );
								    */
								    npve = compare_uv(puv1, data->uv[l], uv_threshold);								  
								    if(  (npve != 0 ) && (wave.eff_wave[k] == data->uv[l].wavelength ) )
								      {
									data->bsref[data->nbis].ab.uvpnt = l;
									data->bsref[data->nbis].ab.sign = 1;
									
									/* conjugated ref if u -ve */
									if(npve == -1)
										{
										  data->bsref[data->nbis].ab.sign = -1;
										}
									uvexists = 1;
									break; /*so that first match is referenced */
								      }
								  }

								if(uvexists == 0)
								  {
								    //								    printf("Warning, orphan bispectrum %d uv1\n", data->nuv);
								  
									/* create new uv point */
									data->uv[data->nuv].u = puv1.u;
									data->uv[data->nuv].v = puv1.v;
									data->uv[data->nuv].wavelength = wave.eff_wave[k];
									data->uv[data->nuv].bandwidth = wave.eff_band[k];
									data->bsref[data->nbis].ab.uvpnt = l;
									data->bsref[data->nbis].ab.sign = 1;

									/* conjugate if u -ve */
									if(data->uv[data->nuv].u<0.0)
									{
										data->uv[data->nuv].u = -puv1.u;
										data->uv[data->nuv].v = -puv1.v;
										data->bsref[data->nbis].ab.sign = -1;
									}
									data->nuv++;
								}

								/*uv2*/
								uvexists = 0;
								for(l=0; l<data->nuv; l++)
								{
								  npve = compare_uv(puv2, data->uv[l], uv_threshold);
								  if(  (npve != 0 ) && (wave.eff_wave[k] == data->uv[l].wavelength ) )
									{
										data->bsref[data->nbis].bc.uvpnt = l;
										data->bsref[data->nbis].bc.sign = 1;

										/* conjugated ref if u -ve */
										if(npve == -1)
										{
											data->bsref[data->nbis].bc.sign = -1;
										}
										uvexists = 1;
										break;
									}
								}

								if(uvexists == 0)
								{
								  //printf("Warning, orphan bispectrum %d uv2\n", data->nuv);
									/* create new uv point */
									data->uv[data->nuv].u = puv2.u;
									data->uv[data->nuv].v = puv2.v;
									data->uv[data->nuv].wavelength = wave.eff_wave[k];
									data->uv[data->nuv].bandwidth = wave.eff_band[k];
									data->bsref[data->nbis].bc.uvpnt = l;
									data->bsref[data->nbis].bc.sign = 1;

									/* conjugate if u -ve */
									if(data->uv[data->nuv].u<0.0)
									{
										data->uv[data->nuv].u = -puv2.u;
										data->uv[data->nuv].v = -puv2.v;
										data->bsref[data->nbis].bc.sign = -1;
									}
									data->nuv++;
								}

								/*uv3 = (-uv2-uv2)*/
								uvexists = 0;
								for(l=0; l<data->nuv; l++)
								{
								  npve = compare_uv(puv3, data->uv[l], uv_threshold);
								  if(  (npve != 0 ) && (wave.eff_wave[k] == data->uv[l].wavelength ) )
									    {
										data->bsref[data->nbis].ca.uvpnt = l;
										data->bsref[data->nbis].ca.sign = 1;

										/* conjugated ref if u -ve */
										if(npve == -1)
										{
											data->bsref[data->nbis].ca.sign = -1;
										}
										uvexists = 1;
										break;
									}
								}

								if(uvexists == 0)
								  {
								    //printf("Warning, orphan bispectrum %d uv3\n", data->nuv);
									/* create new uv point */
									data->uv[data->nuv].u = puv3.u;
									data->uv[data->nuv].v = puv3.v;
									data->uv[data->nuv].wavelength = wave.eff_wave[k];
									data->uv[data->nuv].bandwidth = wave.eff_band[k];
									data->bsref[data->nbis].ca.uvpnt = l;
									data->bsref[data->nbis].ca.sign = 1;

									/* conjugate if u -ve */
									if(data->uv[data->nuv].u<0.0)
									{
										data->uv[data->nuv].u = -puv3.u;
										data->uv[data->nuv].v = -puv3.v;
										data->bsref[data->nbis].ca.sign = -1;
									}
									data->nuv++;
								}

								data->nbis++;
							}
						}
					}
					  printf("\n");
				}
			}
			/* free memory */
			if(*status == 0)
			{
				free_oi_wavelength(&wave);
				free_oi_t3(&t3);
			}
		}
		*status = 0;
	}

	/* ERROR HANDLING */
	return *status;
}

int get_oi_fits_selection(RECONST_PARAMETERS* reconst_parameters, int* status)
{
	/*
	* ROUTINE FOR READING OI_FITS INFO.
	*/
	/* Declare data structure variables*/
	/*     OI-FITS */
	/* oi_array array; */
	oi_target targets;
	oi_wavelength wave;
	oi_vis2 vis2;
	oi_t3 t3;
	/* Data description */
	int nv2tab=0, nt3tab=0;
	int phu;
	/* Declare other variables */

	char comment[FLEN_COMMENT];
	char extname[FLEN_VALUE];
	char zerostring[FLEN_VALUE];
	char commstring[100];
	int hdutype;
	int i,k;
	int nhu=0;
	int tmpi, dummy;
	fitsfile *fptr;


	/* If error do nothing */
	if(*status) return *status;

	/* Initialise */
	reconst_parameters->numins = 0;
	for(k=0; k<FLEN_VALUE-1; k++)zerostring[k] = ' ';
	zerostring[FLEN_VALUE-1]='\0';

	/* Read fits file */
	fits_open_file(&fptr, reconst_parameters->datafile, READONLY, status);
	if(*status) {
		fits_report_error(stderr, *status);
		exit(1);
	}

	/* GET NO OF HEADER UNITS */
	fits_get_num_hdus(fptr,&nhu,status);
	/* PRINT HEADER UNIT LABELS */
	if(*status == 0)
	{
		printf("Reading unit labels:\t");
		for(i = 1; i<(nhu+1); i++)
		{
			fits_movabs_hdu(fptr,i,&hdutype,status);
			if (hdutype == BINARY_TBL) {
				fits_read_key(fptr, TSTRING, "EXTNAME", extname, comment, status);
				printf("%s ",extname);
			}
		}

	}
	/* GET TARGETS AND SELECT */
	if(*status == 0)
	{
		fits_movabs_hdu(fptr,1,NULL,status);
		read_oi_target(fptr,&targets,status);

		for(i=0; i<(targets.ntarget); i++)
		{
			printf("\nTarget id/name:\t\t%d/%s\n",targets.targ[i].target_id,targets.targ[i].target);
		}

		if(targets.ntarget>1)
		{
			AGAIN1:
			printf("\nSELECT AN ID: ");
			dummy = scanf("%d",&reconst_parameters->target_id);
			tmpi = 0;
			for(i=0; i<targets.ntarget; i++)
			{
				if(targets.targ[i].target_id == reconst_parameters->target_id)
				{
					tmpi=1;
					strcpy(reconst_parameters->targetname, targets.targ[i].target);
				}
			}
			if(tmpi==0)goto AGAIN1;
		}
		else
		{
			printf("Auto selecting the only target \"%s\".\n",targets.targ[0].target);
			reconst_parameters->target_id = targets.targ[0].target_id;
			strcpy(reconst_parameters->targetname,zerostring);
			strcpy(reconst_parameters->targetname,targets.targ[0].target);
		}
		printf("\n");
		/* free memory */
		if(*status == 0)
		{
			free_oi_target(&targets);
		}
	}

	/* PRINT AVAILABLE DATA ON SELECTED TARGET */
	if(*status == 0)
	{
		/* V2 TABLES */
		printf("POWERSPECTRUM TABLES\n");
		printf("#\tDate\t\tArray\t\t\tInstrument\t\tNrec/Nwav\n");
		fits_movabs_hdu(fptr,1,NULL,status);
		while(*status == 0)
		{
			read_next_oi_vis2(fptr, &vis2, status);
			fits_get_hdu_num(fptr, &phu);
			if(*status==0)
			{
				if(vis2.record[0].target_id == reconst_parameters->target_id)
				{
					nv2tab++;
					printf("%-6.3d\t%-14.11s\t%-20.20s\t%-20.20s\t%ld/%d\n",nv2tab,vis2.date_obs,vis2.arrname,
					vis2.insname,vis2.numrec,vis2.nwave);
					/* Check if need to register new array for this one */
					tmpi = 0;
					for(i=0; i<reconst_parameters->numins; i++)
					{
						tmpi += ( (int)strcmp(vis2.insname,reconst_parameters->insname[i]) == 0) ;
					}

					if((tmpi == 0)||(reconst_parameters->numins==0))
					{
						strcpy(reconst_parameters->insname[reconst_parameters->numins], zerostring);
						strcpy(reconst_parameters->insname[reconst_parameters->numins], vis2.insname);
						reconst_parameters->numins++;
					}
				}
			}

			/* free memory */
			if(*status == 0)
			{
				free_oi_vis2(&vis2);
			}
		}
		phu = 1;
		*status = 0;
		if(nv2tab<0)printf("\nNo powerspectrum data available for \"%s\"\n",reconst_parameters->targetname);
	}


	/* T3 TABLES */
	if(*status == 0)
	{
		printf("\n");
		printf("BISPECTRUM TABLES\n");
		printf("#\tDate\t\tArray\t\t\tInstrument\t\tNrec/Nwav\n");
		fits_movabs_hdu(fptr,1,NULL,status);

		while(*status == 0)
		{
			read_next_oi_t3(fptr, &t3, status);
			if(*status==0)
			{
				if(t3.record[0].target_id == reconst_parameters->target_id)
				{
					nt3tab++;
					printf("%-6.3d\t%-14.11s\t%-20.20s\t%-20.20s\t%ld/%d\n",nt3tab,
					t3.date_obs,t3.arrname,t3.insname,t3.numrec,t3.nwave);
					/* Check if need to register new array for this one */
					tmpi = 0;
					for(i=0; i<reconst_parameters->numins; i++)
					{
					  tmpi += ( (int)strcmp(t3.insname,reconst_parameters->insname[i]) == 0) ;
					}
					if((tmpi == 0)||(reconst_parameters->numins==0))
					{
						strcpy(reconst_parameters->insname[reconst_parameters->numins],zerostring);
						strcpy(reconst_parameters->insname[reconst_parameters->numins],t3.insname);
						reconst_parameters->numins++;
					}
				}
			}
			if(*status == 0)
			{
				free_oi_t3(&t3);
			}
		}
		phu = 1;
		*status = 0;
		if(nv2tab<0)printf("\nNo bispectrum data available for \"%s\"\n",reconst_parameters->targetname);
	}


	/* CHANNELS IN INSTRUMENTS */
	if(*status == 0)
	  {
	    printf("\nINSTRUMENT SPECTRAL CHANNELS\n");
	    printf("#\tInstrument\t\tChannel_id\tBand/Bandwidth (nm)\n");
	    for(i=0; i<reconst_parameters->numins; i++)
	      {
		fits_movabs_hdu(fptr,1,NULL,status);
		/* Read wave table */
		read_oi_wavelength(fptr, reconst_parameters->insname[i], &wave, status);
		/* Display wave table */
		for(k=0; k<wave.nwave; k++)
		  {
		    if(k==0)printf("%-6d\t%-25.20s",i,reconst_parameters->insname[i]);
		    else printf("%-6.6s\t%-25.20s",zerostring,zerostring);
		    
		    printf("%-3.3d_%-14.3d\t%.0f/%.0f\n",i,k,(wave.eff_wave[k]*billion),(wave.eff_band[k]*billion));
		    
		    /*
		      if (strlen(reconst_parameters->insname[i]) > 0)
		      {
		      read_oi_array(fptr, reconst_parameters->insname[i], &array, status);
		      printf("Note: array %s has %d telescopes\n", vis2.arrname, array.nelement);
		      reconst_parameters->ntelescopes = array.nelement;
		      }
		    */
		  }
		/* free memory */
		if(*status == 0)
		  {
		    free_oi_wavelength(&wave);
		  }
	      }

	  AGAIN2: //waveband selection
	    if( ( reconst_parameters->minband < 0.) || (reconst_parameters->maxband <= reconst_parameters->minband) ) // waveband not set (well) externally
	      {
		if(wave.nwave > 1)
		  {
		    printf("Select a wavelength range (default value = 1 50000) :");
		    if( fgets(commstring,100,stdin) == NULL)
		      {
			printf("Internal error when getting wavelength range\n");
			  getchar();
			  }
		    tmpi = sscanf(commstring,"%f %f", &reconst_parameters->minband, &reconst_parameters->maxband);
		  }
		else
		  {
		    tmpi = -1;
		    printf("Only one spectral channel. ");
		  }
		
		if(tmpi == 2)
		  {
		    if((reconst_parameters->maxband <= reconst_parameters->minband)||(reconst_parameters->minband < 0.0))
		      {
			printf("Invalid band selection!\n");
			goto AGAIN2;
		      }
		  }
		else if(tmpi == -1)
		  {
		    printf("Automatic selection of the full channel\n");
		    reconst_parameters->minband = 1. ; /* (wave.eff_wave[0]-wave.eff_band[0]/2.)*billion ; */
		    reconst_parameters->maxband = 50000. ; /*(wave.eff_wave[0]+wave.eff_band[0]/2.)*billion ; */
		  }
		else
		  {
		    printf("Invalid band selection!\n");
		    goto AGAIN2;
		  }
	      }
	    
	  }
	
	/* Count number of vis2 and t3 available in this range for this target */
	if(*status==0)
	  {
	    /* powerspectrum */
	    reconst_parameters->numvis2 = 0;
	    fits_movabs_hdu(fptr,1,NULL,status);
	    while(*status==0)
	      {
		read_next_oi_vis2(fptr, &vis2, status);
		fits_get_hdu_num(fptr, &phu);
		read_oi_wavelength(fptr, vis2.insname, &wave, status);
		fits_movabs_hdu(fptr,phu,NULL,status);
		
		if(*status==0)
		  {
		    if(vis2.record[0].target_id == reconst_parameters->target_id)
		      {
			for(i=0; i<vis2.numrec; i++)
			  {
			    for(k=0; k<vis2.nwave; k++)
			      {
				if(((wave.eff_wave[k]*billion)>reconst_parameters->minband)&&((wave.eff_wave[k]*billion)<reconst_parameters->maxband)&&(!(vis2.record[i].flag[k])))
				  {
				    reconst_parameters->numvis2++;
				  }
			      }
			  }
		      }
		  }
		/* free memory */
		if(*status == 0)
		  {
		    free_oi_vis2(&vis2);
		    free_oi_wavelength(&wave);
		  }
	      }
	    *status = 0;
	    
	    /* bispectrum */
	    reconst_parameters->numt3 = 0;
	    fits_movabs_hdu(fptr,1,NULL,status);
	    while(*status==0)
	      {
		read_next_oi_t3(fptr, &t3, status);
		fits_get_hdu_num(fptr, &phu);
		read_oi_wavelength(fptr, t3.insname, &wave, status);
		fits_movabs_hdu(fptr,phu,NULL,status);
		
		if(*status==0)
			{
			  if(t3.record[0].target_id == reconst_parameters->target_id)
			    {
			      for(i=0; i<t3.numrec; i++)
				{
				  for(k=0; k<t3.nwave; k++)
				    {
				      if( ((wave.eff_wave[k]*billion)>reconst_parameters->minband)
					  &&((wave.eff_wave[k]*billion)<reconst_parameters->maxband)
					  && (!(t3.record[i].flag[k]))  )
					{
					  reconst_parameters->numt3++;
							}
				    }
				}
			    }
			}
		/* free memory */
		if(*status == 0)
		  {
				free_oi_t3(&t3);
				free_oi_wavelength(&wave);
		  }
	      }
		*status = 0;
		if((reconst_parameters->numt3==0)&&(reconst_parameters->numvis2==0))
		{
		  printf("Error: no data available within the selected waveband limits\n");
		  reconst_parameters->minband= -1.;
		  reconst_parameters->maxband= -1.;
		  getchar();
		  goto AGAIN2;
		}
		printf("Found %ld powerspectrum and %ld bispectrum points between %.0f and %.0f nm.\n\n",
						reconst_parameters->numvis2,reconst_parameters->numt3, reconst_parameters->minband, reconst_parameters->maxband);
		
	  }
	
	/* CLOSE FILE */
	fits_close_file(fptr, status);	
	/* ERROR HANDLING */
	//	printf("DEBUG\n");
	return *status;
}

void free_oi_target(oi_target *targets)
{
	free(targets->targ);
}

void free_oi_wavelength(oi_wavelength *wave)
{
	free(wave->eff_wave);
	free(wave->eff_band);
}

void free_oi_vis2(oi_vis2 *vis2)
{
	int i;

	for(i=0; i<vis2->numrec; i++)
	{
		free(vis2->record[i].vis2data);
		free(vis2->record[i].vis2err);
		free(vis2->record[i].flag);
	}
	free(vis2->record);
}

void free_oi_t3(oi_t3 *t3)
{
	int i;

	for(i=0; i<t3->numrec; i++)
	{
		free(t3->record[i].t3amp);
		free(t3->record[i].t3amperr);
		free(t3->record[i].t3phi);
		free(t3->record[i].t3phierr);
		free(t3->record[i].flag);
	}
	free(t3->record);
}

void free_oi_data(oi_data *data)
{
	free(data->pow);
	free(data->powerr);
	free(data->bisamp);
	free(data->bisamperr);
	free(data->bisphs);
	free(data->bisphserr);
	free(data->uv);
	free(data->bsref);
}

int count_redundant_bsuv(oi_bsref *bsref, int nbs)
{
	int i,k;
	int c = 0;
	int fab,fbc,fca;

	for(i=0; i<nbs; i++)
	{
		fab = 0;
		fbc = 0;
		fca = 0;
		for(k=0; k<nbs; k++)
		{
			if( (bsref[i].ab.uvpnt == bsref[k].ab.uvpnt) ||
				(bsref[i].ab.uvpnt == bsref[k].bc.uvpnt) ||
			(bsref[i].ab.uvpnt == bsref[k].ca.uvpnt)    )fab++;

			if( (bsref[i].bc.uvpnt == bsref[k].ab.uvpnt) ||
				(bsref[i].bc.uvpnt == bsref[k].bc.uvpnt) ||
			(bsref[i].bc.uvpnt == bsref[k].ca.uvpnt)    )fbc++;

			if( (bsref[i].ca.uvpnt == bsref[k].ab.uvpnt) ||
				(bsref[i].ca.uvpnt == bsref[k].bc.uvpnt) ||
			(bsref[i].ca.uvpnt == bsref[k].ca.uvpnt)    )fca++;
		}

		if(fab>1)c++;
		if(fbc>1)c++;
		if(fca>1)c++;
	}

	return c;
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
