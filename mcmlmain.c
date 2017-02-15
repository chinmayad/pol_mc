/***********************************************************
** National Institute of technology Karnataka
** Indian Institute of Science, Bangalore, Karnataka, India 560012.
** 2014.
*
*	Launch, move, and record photon weight.
*
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al and Dr.Manojit et al.
**  This code is modified by Chinmaya,Rajat hebbar,Rakesh M,Sharath,Rajan K. 

 *
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in
 *	Think C. If 1, remember to turn on "Generate profiler
 *	calls" in the options menu.
 ****/
#define THINKCPROFILER 0



/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif
#define	NN      	100
#include "mcml.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "array.h"
#include "complex.h"
#include "mie.h"
#include "nrutil.h"
#include <time.h>


/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *);
void FreeData(InputStruct, OutStruct *);
double Rspecular(LayerStruct * );
void LaunchPhoton(double, LayerStruct *, PhotonStruct *,Stokes_struct *);
void LaunchPhotonObj(double, InputStruct *, PhotonStruct *,Stokes_struct *, long);	/**	Included by Vijitha et al **/
void HopDropSpin(InputStruct  *,PhotonStruct *,Stokes_struct *,OutStruct *);
void SumScaleResult(InputStruct, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *);
void 	rotSphi(double* S, double phi, double* S2);
double	RandomGen(char Type, long Seed, long *Status);
void FileWrite(OutStruct * Out_Ptr);
void file(PhotonStruct * Photon_Ptr,Stokes_struct * Stokes_Ptr,long int n);

/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on
 *	screen, return the real time since F=0.
 *
 *	If F = 2, same as F=1 except no printing.
 *
 *	Note that clock() and time() return user time and real
 *	time respectively.
 *	User time is whatever the system allocates to the
 *	running of the program;
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *
 *	clock() only hold 16 bit integer, which is about 32768
 *	clock ticks.
 ****/
time_t PunchTime(char F, char *Msg)
{
#if GNUCC
  return(0);
#else
  static clock_t ut0;	/* user time reference. */
  static time_t  rt0;	/* real time reference. */
  double secs;
  char s[STRLEN];

  if(F==0) {
    ut0 = clock();
    rt0 = time(NULL);
    return(0);
  }
  else if(F==1)  {
    secs = (clock() - ut0)/(double)CLOCKS_PER_SEC;
    if (secs<0) secs=0;	/* clock() can overflow. */
    sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n",
	    secs, secs/3600.0, Msg);
    puts(s);
    strcpy(Msg, s);
    return(difftime(time(NULL), rt0));
  }
  else if(F==2) return(difftime(time(NULL), rt0));
  else return(0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)
{
  time_t now, done_time;
  struct tm *date;
  char s[80];

  now = time(NULL);
  date = localtime(&now);
  strftime(s, 80, "%H:%M %x", date);
  printf("Now %s, ", s);

  done_time = now +
			(time_t) (PunchTime(2,"")/(double)P1*(Pt-P1));
  date = localtime(&done_time);
  strftime(s, 80, "%H:%M %x", date);
  printf("End %s\n", s);
}

/***********************************************************
 *	Report time and write results.
 ****/
void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
  char time_report[STRLEN];

  strcpy(time_report, " Simulation time of this run.");
//  PunchTime(1, time_report);

 // SumScaleResult(In_Parm, &Out_Parm);
 // WriteResult(In_Parm, Out_Parm, time_report);
}

/***********************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc,
					  char * argv[],
					  char * input_filename)
{
  if(argc>=2) {			/* filename in command line */
    strcpy(input_filename, argv[1]);
  }
  else
    input_filename[0] = '\0';
}


/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(short NumRuns, InputStruct *In_Ptr,Stokes_struct * Stokes_Ptr)
{

  register long i_photon;
	/* index to photon. register for speed.*/
  OutStruct out_parm;		/* distribution of photons.*/
  PhotonStruct photon;
  long num_photons = In_Ptr->num_photons, photon_rep=50000; /** change this number to change the display of the number of photons left for run **/
  long temp = num_photons;
  //photon.pathlength=(PhotonStruct *)malloc(sizeof(double)*num_photons);
  photon.pathlength=(double *)malloc(sizeof(double)*num_photons);
    photon.pathlengthinobj=(double *)malloc(sizeof(double)*num_photons);

  short i,j;



  for(i=0;i<200;i++)
  {
	  for(j=0;j<200;j++)
	  {
		  out_parm.A_xz[i][j]=0;
	  }
  }

 // for(i=0;i<num_photons;i++)
//	  photon.pathlength[i]=0;

#if THINKCPROFILER
  InitProfile(200,200); cecho2file("prof.rpt",0, stdout);
#endif

  InitOutputData(*In_Ptr, &out_parm);
  out_parm.Rsp = Rspecular(In_Ptr->layerspecs);
  i_photon = num_photons;
  photon.i_photon = num_photons;
  photon.n_photon = num_photons;
  PunchTime(0, "");
  for(i=0;i<200;i++)
    {
  	  for(j=0;j<200;j++)
  	  {

  		  out_parm.Rd_xy[i][j]=0;
  		  out_parm.Tt_xy[i][j]=0;
  		  out_parm.IR[i][j]=0;
  		  out_parm.QR[i][j]=0;
  		  out_parm.UR[i][j]=0;
  		  out_parm.VR[i][j]=0;
  		  out_parm.ITT[i][j]=0;
  		  out_parm.QTT[i][j]=0;
  		  out_parm.UTT[i][j]=0;
  		  out_parm.VTT[i][j]=0;
  	  }
  	  out_parm.P_z[i]=0;
    }

  out_parm.IR_1=0;
  out_parm.QR_1=0;
  out_parm.UR_1=0;
  out_parm.VR_1=0;


  do {
    if(temp - i_photon == photon_rep) {
	  temp = i_photon;
      printf("%ld photons & %hd runs left, ", i_photon, NumRuns);
      PredictDoneTime(num_photons - i_photon, num_photons);
     // photon_rep *= 10;
    }
	if (In_Ptr->objLayer == 0 || In_Ptr->objCode == 0) {		/** launch photon when there is no object **/
		LaunchPhoton(out_parm.Rsp, In_Ptr->layerspecs, &photon,Stokes_Ptr);
	}
	else {								/** launch photon when there is object with inobj pointer **/
		LaunchPhotonObj(out_parm.Rsp, In_Ptr, &photon,Stokes_Ptr, i_photon);
	}
	do  HopDropSpin(In_Ptr, &photon, Stokes_Ptr, &out_parm);
    while (!photon.dead);

	printf("%ld\n",i_photon);
	--photon.i_photon;
  } while(--i_photon);

  printf("\n IR_1 = %f\n QR_1 = %f\n UR_1 = %f\n VR_1 = %f \n",out_parm.IR_1,out_parm.QR_1,out_parm.UR_1,out_parm.VR_1);

#if THINKCPROFILER
  exit(0);
#endif

  FileWrite(&out_parm);
  file(&photon,Stokes_Ptr,num_photons);

  // ReportResult(*In_Ptr, out_parm);
  //FreeData(*In_Ptr, &out_parm);
}

/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/
int
main(int argc, char *argv[])
{
	double pi = 3.1415926535897932384;

			/* Mie theory stuff */
			double radius[10],lambda, A[10];
			long nangles,i;
			struct complex m;
			struct complex*s1=NULL;
			struct complex*s2=NULL;
			Stokes_struct * Stokes_Ptr=(Stokes_struct *)malloc(sizeof(Stokes_struct));
			double *mu=NULL;
			double x[10],qext,qsca,qback,g, rho[10], vol[10],dy, dx, hw;
			double nre_p[10], nim_p, nre_med[10], nim_med;
			FILE *target;
			double jjj;
			double mua,Nphotons;
			double mus[10];
			/* E field stuff */
			double	phi, theta,I,I0;
			int 	ithedeg;
			double  me[10];


			int     MM;
				double  W,absorb;  /* photon weight */
				double  slabsize;
				int 	j,ix,iy;
				double  cos22,sin22,costheta,sini,cosi;


				//double	*S;     	/* */
				//	double	*S0;     	/* */
				//	double	*S2;     	/* */
					/*double	*s11=NULL;
					double	*s12=NULL;
					double	*s33=NULL;
					double	*s43=NULL;*/


					MM = NN - 1;
	/*
					double S[4];
					double S0[4];
					double S2[4];

					double IR[99][99];
					double QR[99][99];
					double UR[99][99];
					double VR[99][99];

					double mu[1000];
					double s11[1000];
					double s12[1000];
					double s33[1000];
					double s43[1000];
	*/
					//S      = new_darray(4);
					//S0     = new_darray(4);
					//S2     = new_darray(4);/* dummy S*/

					//IQUV   = new_darray(4);

					/* CHOOSE MIE SCATTERING parameters */
					//radius  	= 2.03/2; /* microns */
					lambda 		= 0.6328; /* microns */
					//rho 		= 1.152e-4;/*Dilution 1*/
					Nphotons	= 1e6;
					mua 		= 0.0; /*�a  */

					/* ------------------------*/
					//nre_p   	= 1.59;
					nim_p  	 	= 0;
					//nre_med		= 1.33;
					nim_med 	= 0.0;
					nangles 	= 1000;


					mu  = new_darray(nangles);
					s1  = new_carray(nangles);
					s2  = new_carray(nangles);
				/*	s11 = new_darray(nangles);
					s12 = new_darray(nangles);
					s33 = new_darray(nangles);
					s43 = new_darray(nangles);*/

					nre_p[1]   	= 1.36;
					nre_p[2]   	= 1.7;
					nre_p[3]   	= 1.43;
					nre_p[4]   	= 1.41;
					nre_p[5]   	= 1.36;
					nre_p[6]   	= 1.41;
					nre_p[7]   	= 1.39;

					nre_med[1]		= 1.5;
					nre_med[2]		= 1.34;
					nre_med[3]		= 1.4;
					nre_med[4]		= 1.39;
					nre_med[5]		= 1.4;
					nre_med[6]		= 1.38;
					nre_med[7]		= 1.44;


					radius[1]  	= 0.53/2;/*microns*/
					radius[2]  	= 0.471/2;
					radius[3]  	= 0.68/2;
					radius[4]  	= 1.04/2;
					radius[5]  	= 0.44/2;
					radius[6]  	= 1.04/2;
					radius[7]  	= 0.37/2;

					rho[1]=2.2;
					rho[2]=0.2;
					rho[3]=4.4;
					rho[4]=2;
					rho[5]=13;
					rho[6]=0.76;
					rho[7]=3.5;

					for(j=1;j<8;j++){
						x[j]    = 2*pi*radius[j]/(lambda/nre_med[j]);
						vol[j]  = 4.0/3*pi*radius[j]*radius[j]*radius[j];
						A[j]    = pi*radius[j]*radius[j];
						//me[j]	= nre_p[j]/nre_med[j];
					}


					//m.re = nre_p/nre_med;
					//m.im = 0.0;

					for(i=0;i<=nangles;i++)
						mu[i] = cos(pi*i/nangles);
						/*s11=new_darray(nangles);
						s12=new_darray(nangles);
						s33=new_darray(nangles);
						s43=new_darray(nangles);*/
						s1=new_carray(nangles);
						s2=new_carray(nangles);

					//	mus 	= qsca*A*rho*1e4; /* Mus is in cm^-1 */
						//musp 	= mus*(1-g);/* [cm^-1] */
						//albedo 	= mus/(mus + mua);
						//free_darray(mu);
						//Mie(x,m,mu,nangles,s1,s2,&qext,&qsca,&qback,&g); /* <---- Call Mie program ----- */

						/*for(j=0;j<10;j++){
						for(i=0;i<=nangles;++i){
							Stokes_Ptr->s11[j][i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) + 0.5*cabbs(s1[i])*cabbs(s1[i]);
															Stokes_Ptr->s12[j][i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) - 0.5*cabbs(s1[i])*cabbs(s1[i]);
																			Stokes_Ptr->s33[j][i] = (cmul(conj(s1[i]),s2[i])).re;
																			Stokes_Ptr->s43[j][i] = (cmul(conj(s1[i]),s2[i])).im;
						}
						}*/

						for(j=1;j<8;j++)
												{
												m.re=nre_p[j]/nre_med[j];
												m.im=0;
												Mie(x[j],m,mu,nangles,s1,s2,&qext,&qsca,&qback,&g); /* <---- Call Mie program ----- */

												for(i=0;i<=nangles;++i){
													Stokes_Ptr->s11[j][i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) + 0.5*cabbs(s1[i])*cabbs(s1[i]);
													Stokes_Ptr->s12[j][i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) - 0.5*cabbs(s1[i])*cabbs(s1[i]);
													Stokes_Ptr->s33[j][i] = (cmul(conj(s1[i]),s2[i])).re;
													Stokes_Ptr->s43[j][i] = (cmul(conj(s1[i]),s2[i])).im;
												/*printf("%5.5f\t %5.5f\t %5.5f\t %5.5f\n",s11[i],s12[i],s33[i],s43[i]);*/

													}
												mus[j] 	= qsca*A[j]*rho[j]*1e4;
							}

							//hw 			= 7/mus; /* [cm] , maximum range in x and y for output. */
						//	dx 			= 2.0*hw/NN;
							//dy 			= 2.0*hw/NN;

						/******** MONTE CARLO *******/
  char input_filename[STRLEN];
  FILE *input_file_ptr;
  short num_runs;	/* number of independent runs. */
  InputStruct in_parm;

  ShowVersion("Version 1.0.0, 2014");	/** Updated by Vijitha et al **/
  GetFnameFromArgv(argc, argv, input_filename);
  input_file_ptr = GetFile(input_filename);
  CheckParm(input_file_ptr, &in_parm);
  num_runs = ReadNumRuns(input_file_ptr);

  while(num_runs--)  {
    ReadParm(input_file_ptr, &in_parm);
	DoOneRun(num_runs, &in_parm,Stokes_Ptr);
  }
  fclose(input_file_ptr);
  return(0);
}

void FileWrite(OutStruct * Out_Ptr)
{
	short i,j;
	FILE * file;
	  file=fopen("ir.dat","w");
	  for(i=0;i<200;i++)
	  {
		  for(j=0;j<200;j++)
		  {
			  fprintf(file,"%f ",Out_Ptr->IR[i][j]);
		  }
		  fprintf(file,"\n");
	  }
	  fclose(file);

	  FILE * file1;
	   file1=fopen("qr.dat","w");
	   for(i=0;i<200;i++)
	   {
	 	  for(j=0;j<200;j++)
	 	  {
	 		  fprintf(file1,"%f ",Out_Ptr->QR[i][j]);
	 	  }
	 	  fprintf(file1,"\n");
	   }
	   fclose(file1);

	   FILE * file2;
	      file2=fopen("ur.dat","w");
	      for(i=0;i<200;i++)
	      {
	    	  for(j=0;j<200;j++)
	    	  {
	    		  fprintf(file2,"%f ",Out_Ptr->UR[i][j]);
	    	  }
	    	  fprintf(file2,"\n");
	      }
	      fclose(file2);

	      FILE * file3;
	      	      file3=fopen("vr.dat","w");
	      	      for(i=0;i<200;i++)
	      	      {
	      	    	  for(j=0;j<200;j++)
	      	    	  {
	      	    		  fprintf(file3,"%f ",Out_Ptr->VR[i][j]);
	      	    	  }
	      	    	  fprintf(file3,"\n");
	      	      }
	      	      fclose(file3);


	      FILE * file4;
	      	     file4=fopen("itt.dat","w");
	             for(i=0;i<200;i++)
	      	     {
	      	     	   for(j=0;j<200;j++)
	      	      	   {
	      	    	 	   fprintf(file4,"%f ",Out_Ptr->ITT[i][j]);
	      	    	   }
	      	    	   fprintf(file4,"\n");
	      	     }
	      	     fclose(file4);

	      FILE * file5;
	   	      file5=fopen("qtt.dat","w");
	   	      for(i=0;i<200;i++)
	   	      {
	   	    	  for(j=0;j<200;j++)
	   	    	  {
	   	    		  fprintf(file5,"%f ",Out_Ptr->QTT[i][j]);
	   	    	  }
	   	    	  fprintf(file5,"\n");
	   	      }
	   	      fclose(file5);

	   	   FILE * file6;
	   	   	      file6=fopen("utt.dat","w");
	   	   	      for(i=0;i<200;i++)
	   	   	      {
	   	   	    	  for(j=0;j<200;j++)
	   	   	    	  {
	   	   	    		  fprintf(file6,"%f ",Out_Ptr->UTT[i][j]);
	   	   	    	  }
	   	   	    	  fprintf(file6,"\n");
	   	   	      }
	   	   	      fclose(file6);
	   	   	FILE * file7;
	   	   		      file7=fopen("vtt.dat","w");
	   	   		      for(i=0;i<200;i++)
	   	   		      {
	   	   		    	  for(j=0;j<200;j++)
	   	   		    	  {
	   	   		    		  fprintf(file7,"%f ",Out_Ptr->VTT[i][j]);
	   	   		    	  }
	   	   		    	  fprintf(file7,"\n");
	   	   		      }
	   	   		      fclose(file7);

	   	   		FILE * file8;
	   	   			   	   		      file8=fopen("P_z.dat","w");

	   	   			   	   		    	  for(j=0;j<200;j++)
	   	   			   	   		    	  {
	   	   			   	   		    		  fprintf(file8,"%f ",Out_Ptr->P_z[j]);
	   	   			   	   		    	  }


	   	   			   	   		      fclose(file8);


	   	FILE * file9;
	   	   	 file9=fopen("xyz.dat","w");

	   	   	 for(i=0;i<200;i++)
	   	   	 {
	   	   		 for(j=0;j<200;j++)
	   	   		 {
	   	   			 fprintf(file9,"%lf ",Out_Ptr->A_xz[i][j]);
	   	   		 }
	   	   		 fprintf(file9,"\n");
	   	   	 }
	   	   	 fclose(file9);


}


void file(PhotonStruct * Photon_Ptr,Stokes_struct * Stokes_Ptr, long int n)
{

	short i,j;
	FILE * file;
	file=fopen("Pathlength.dat","w");
	for(i=0;i<n;i++)
	{
		fprintf(file,"%f\n",Photon_Ptr->pathlength[i]);
	}
	fclose(file);

	FILE * file1;
	file1=fopen("PathlengthInObj.dat","w");
	for(i=0;i<n;i++)
	{
		fprintf(file1,"%f\n",Photon_Ptr->pathlengthinobj[i]);
	}
	fclose(file1);

	FILE * file2;
	file2=fopen("Final_Coordinates.dat","w");
	for(i=0;i<n;i++)
	{
		fprintf(file2,"%f, %f, %f\n",Photon_Ptr->x,Photon_Ptr->y,Photon_Ptr->z);
	}
	fclose(file2);

	FILE * file3;
	file3=fopen("Final_Stokes_Vector.dat","w");
	for(i=0;i<n;i++)
	{
		fprintf(file3,"%f, %f, %f, %f\n",Stokes_Ptr->S[0],Stokes_Ptr->S[1],Stokes_Ptr->S[2],Stokes_Ptr->S[3]);
	}
	fclose(file3);
}
