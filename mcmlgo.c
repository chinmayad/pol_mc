/***********************************************************
** 
** National Institute of technology Karnataka
** Indian Institute of Science, Bangalore, Karnataka, India 560012.
** 2014.
*
*	Launch, move, and record photon weight.
*
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al and Dr.Manojit et al.
**  This code is modified by Chinmaya,Rajat hebbar,Rakesh M,Sharath,Rajan K. 

****/

#include "mcml.h"

#define RandomNum1 (double) RandomGen(1, 0, NULL)

#define InitRandomGen (double) RandomGen(0, 1, NULL)

#define STANDARDTEST 0
/* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0
/* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)
/* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6
/* cosine of about 1.57 - 1e-6 rad. */


//#define InitRandomGen (double) RandomGen(0, 1, NULL)

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

void rotSphi(double* S, double phi, double* S2);
/********************************************************************
*	Compute the Stokes Vector after reflection from a boundary.
*
*	Stokes Vector S2 is reflected and stored in S (since result of rotSphi is in S2).
*
*	****/
void StokesReflect(Stokes_struct * Stokes_Ptr,double i,double * ta)
{
	double ttm,stp,ctp,ctm;
	double S[4],S2[4];
	short j;
	double r=acos(*ta);

		for (j=0;j<4;j++)
			{
				S[j]=Stokes_Ptr->S2[j];
			}


	ttm=tan(i-r);
	stp=sin(i+r);
	ctp=cos(i+r);
	ctm=cos(i-r);

	S2[0]=0.5*ttm*ttm/(stp*stp)*((ctm*ctm+ctp*ctp)*S[0]+(ctm*ctm-ctp*ctp)*S[1]);
	S2[1]=0.5*ttm*ttm/(stp*stp)*((ctm*ctm-ctp*ctp)*S[0]+(ctm*ctm+ctp*ctp)*S[1]);
	S2[2]=0.5*ttm*ttm/(stp*stp)*((-2)*ctp*ctm)*S[2];
	S2[3]=0.5*ttm*ttm/(stp*stp)*((-2)*ctp*ctm)*S[3];


	for (j=0;j<4;j++)
			{
				Stokes_Ptr->S[j]=S2[j]/S2[0];
			}

	if(i==0)											//For normal incidence.
		*ta=1;

}
/********************************************************************
*	Compute the Stokes Vector after transmission through a boundary.
*
*	Stokes Vector S2 is transmitted and stored in S (since result of rotSphi is in S2).
*
*	****/
void StokesTransmit(Stokes_struct * Stokes_Ptr,double i,double * ta,double n1,double n2)
{
	double k;
	double stp,ctm,s2i,s2r;
	double S[4],S2[4];
	short j;
	double r = acos(*ta);

	if(i==0)											//For normal incidence.
	{
		for (j=0;j<4;j++)
					{
						Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
					}
		*ta=1;
	}
	else												//For regular incidence.
	{
		for (j=0;j<4;j++)
			{
				S[j]=Stokes_Ptr->S2[j];
			}

	stp=sin(i+r);
	ctm=cos(i-r);
	s2i=sin(2*i);
	s2r=sin(2*r);
	k=(n2/n1)*s2r*s2i/(2*stp*stp*ctm*ctm);

	S2[0]=k*((ctm*ctm+1)*S[0]+(ctm*ctm-1)*S[1]);
	S2[1]=k*((ctm*ctm-1)*S[0]+(ctm*ctm+1)*S[1]);
	S2[2]=k*2*ctm*S[2];
	S2[3]=k*2*ctm*S[3];
	for (j=0;j<4;j++)
			{
				Stokes_Ptr->S[j]=S2[j]/S2[0];
			}
	}
}

double
RandomGen(char Type, long Seed, long *Status){
  static long i1, i2, ma[56];   /* ma[0] is not used. */
  long        mj, mk;
  short       i, ii;

  if (Type == 0) {              /* set seed. */
    mj = MSEED - (Seed < 0 ? -Seed : Seed);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (ii = 1; ii <= 4; ii++)
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    i1 = 0;
    i2 = 31;
  } else if (Type == 1) {       /* get a number. */
    if (++i1 == 56)
      i1 = 1;
    if (++i2 == 56)
      i2 = 1;
    mj = ma[i1] - ma[i2];
    if (mj < MZ)
      mj += MBIG;
    ma[i1] = mj;
    return (mj * FAC);
  } else if (Type == 2) {       /* get status. */
    for (i = 0; i < 55; i++)
      Status[i] = ma[i + 1];
    Status[55] = i1;
    Status[56] = i2;
  } else if (Type == 3) {       /* restore status. */
    for (i = 0; i < 55; i++)
      ma[i + 1] = Status[i];
    i1 = Status[55];
    i2 = Status[56];
  } else
    puts("Wrong parameter to RandomGen().");
  return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/***********************************************************
*	A random number generator from Numerical Recipes in C.
****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

float ran3(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
			inext=0;
			inextp=31;
			*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return (double)mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
*	Generate a random number between 0 and 1.  Take a
*	number as seed the first time entering the function.
*	The seed is limited to 1<<15.
*	We found that when idum is too large, ran3 may return
*	numbers beyond 0 and 1.
****/
double RandomNum(void)
{
	static Boolean first_time=1;
	static int idum;	/* seed for ran3. */

	if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
		idum = - 1;
#else
		idum = -(int)time(NULL)%(1<<15);
		/* use 16-bit integer as the seed. */
#endif
		ran3(&idum);
		first_time = 0;
		idum = 1;
	}

	return( (double)ran3(&idum) );
}

/***********************************************************
*	Compute the specular reflection.
*
*	If the first layer is a turbid medium, use the Fresnel
*	reflection from the boundary of the first layer as the
*	specular reflectance.
*
*	If the first layer is glass, multiple reflections in
*	the first layer is considered to get the specular
*	reflectance.
*
*	The subroutine assumes the Layerspecs array is correctly
*	initialized.
****/
double Rspecular(LayerStruct * Layerspecs_Ptr)
{
	double r1, r2;
	/* direct reflections from the 1st and 2nd layers. */
	double temp;

	temp =(Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
		/(Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
	r1 = temp*temp;

	if((Layerspecs_Ptr[1].mua == 0.0)
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n)
				/(Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);
			r2 = temp*temp;
			r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
	}

	return (r1);
}

/***********************************************************
*	Compute the Fresnel reflectance.
*
*	Make sure that the cosine of the incident angle a1
*	is positive, and the case when the angle is greater
*	than the critical angle is ruled out.
*
* 	Avoid trigonometric function operations as much as
*	possible, because they are computation-intensive.
****/
double RFresnel(Stokes_struct * Stokes_Ptr,double n1,	/* incident refractive index.*/
	double n2,	/* transmit refractive index.*/
	double ca1,	/* cosine of the incident */
	/* angle. 0<a1<90 degrees. */
	double * ca2_Ptr)  /* pointer to the */
	/* cosine of the transmission */
	/* angle. a2>0. */
{
	double r,rs,rp;
	double P,Q;
	Q=Stokes_Ptr->S[1];

	if(n1==n2) {			  	/** matched boundary. **/
		*ca2_Ptr = ca1;
		r = 0.0;
	}
	else if(ca1>COSZERO) {	/** normal incident. **/
		*ca2_Ptr = ca1-COS90D;   		//So that incident and refracted angles are not exactly equal, will affect StokesReflect and StokesTransmit.
		r = (n2-n1)/(n2+n1);
		r *= r;
	}
	else if(ca1<COS90D)  {	/** very slant. **/
		*ca2_Ptr = 0.0;
		r = 1.0;
	}
	else  {			  		/** general. **/
		double sa1, sa2;
		/* sine of the incident and transmission angles. */
		double ca2;

		sa1 = sqrt(1-ca1*ca1);
		sa2 = n1*sa1/n2;
		if(sa2>=1.0) {
			/* double check for total internal reflection. */
			*ca2_Ptr = 0.0;
			r = 1.0;
		}
		else  {
			double cap, cam;	/* cosines of the sum ap or */
			/* difference am of the two */
			/* angles. ap = a1+a2 */
			/* am = a1 - a2. */
			double sap, sam;	/* sines. */
						double tap, tam;    /* tangents */

						*ca2_Ptr = ca2 = sqrt(1-sa2*sa2);

						cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
						cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
						sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
						sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
						tap = sap/cap;
						tam = sam/cam;

						rs = sam*sam/(sap*sap);
						rp = tam*tam/(tap*tap);

						double P=(Q+1)/2;
						r=((1-P)*rs+P*rp);
			/* rearranged for speed. */
		}
	}
	return(r);
}

/***********************************************************
*	Initialize a photon packet.
****/
void LaunchPhoton(double Rspecular,
	LayerStruct  * Layerspecs_Ptr,
	PhotonStruct * Photon_Ptr,
	Stokes_struct * Stokes_Ptr)
{
	short i;
	double a;
	double S0[4]={1,1,0,0};
	Photon_Ptr->w	 	= 1.0 - Rspecular;//removed temporarily, required for MCML
	Photon_Ptr->dead 	= 0;
	Photon_Ptr->layer = 1;
	Photon_Ptr->s	= 0;
	Photon_Ptr->sleft= 0;

	Photon_Ptr->x 	= 0.0;
	Photon_Ptr->y	= 0.0;
	Photon_Ptr->z	= 0.0;
	Photon_Ptr->ux	= 0.0;
	Photon_Ptr->uy	= 0.0;
	Photon_Ptr->uz	= 1;
	Photon_Ptr->inObj	= 0;
	Photon_Ptr->pathlength[Photon_Ptr->n_photon-Photon_Ptr->i_photon]=0;

	if(Photon_Ptr->i_photon==Photon_Ptr->n_photon)
	{	InitRandomGen;
	//a=RandomNum1;
	//Photon_Ptr->count++;
	}


	if((Layerspecs_Ptr[1].mua == 0.0)
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			Photon_Ptr->layer 	= 2;
			Photon_Ptr->z	= Layerspecs_Ptr[2].z0;
	}
	for(i=0;i<4;i++)
	    						{
	    							Stokes_Ptr->S[i]=S0[i];
	    						     Stokes_Ptr->S2[i]=0;
	    						}

}


/***********************************************************
**	Initialize a photon packet when there is a object.
****/
void LaunchPhotonObj(double Rspecular,
					 InputStruct * In_Ptr,
					 PhotonStruct * Photon_Ptr,
					 Stokes_struct * Stokes_Ptr,
					 long i_photon)
{
	short i,count=0;
	double a;
	double S0[4]={1,1,0,0};
	double uz1, ni, nt, uz, r;
	short launchConfig;
	long numOfPhotons;
	ni = In_Ptr->layerspecs[0].n;
	nt = In_Ptr->layerspecs[1].n;
	numOfPhotons = In_Ptr->num_photons;

	Photon_Ptr->inObj   = 0;  /** Initialize that the photon is outside object. **/
	Photon_Ptr->w	 	= 1.0 - Rspecular;   //removed temporarily, required for MCML
	Photon_Ptr->dead 	= 0;
	Photon_Ptr->layer   = 1;
	Photon_Ptr->s	    = 0;
	Photon_Ptr->sleft   = 0;

	Photon_Ptr->x 	= 0.0;
	Photon_Ptr->y	= 0.0;
	Photon_Ptr->z	= 0.0;
	Photon_Ptr->ux	= 0.0;
	Photon_Ptr->uy	= 0.0;
	Photon_Ptr->uz	= 1;

	Photon_Ptr->pathlength[Photon_Ptr->n_photon-Photon_Ptr->i_photon]=0;

	if(Photon_Ptr->i_photon==Photon_Ptr->n_photon)
	{

		InitRandomGen ;
		//a=RandomNum1;
		//Photon_Ptr->count++;
	}

	if((In_Ptr->layerspecs[1].mua == 0.0)
		&& (In_Ptr->layerspecs[1].mus == 0.0))  { /* glass layer. */
			Photon_Ptr->layer = 2;
			Photon_Ptr->z = In_Ptr->layerspecs[2].z0;
	}
	for(i=0;i<4;i++)
		    						{
		    							Stokes_Ptr->S[i]=S0[i];
		    						     Stokes_Ptr->S2[i]=0;
		    						}
}

/***********************************************************
*	Choose a new direction for photon propagation by
*	sampling the polar deflection angle theta and the
*	azimuthal angle psi.
*
*	Note:
*  	theta: 0 - pi so sin(theta) is always positive
*  	feel free to use sqrt() for cos(theta).
*
*  	psi:   0 - 2pi
*  	for 0-pi  sin(psi) is +
*  	for pi-2pi sin(psi) is -
*
*  	Rejection is used to choose theta and psi.
*
*  	Stokes vector is updated as in stok1.
****/
void Spin(double g,
	PhotonStruct * Photon_Ptr,
	Stokes_struct * Stokes_Ptr)
{
	double cost, sint;	/* cosine and sine of the */
	/* polar deflection angle theta. */
	double cosp, sinp;	/* cosine and sine of the */
	/* azimuthal angle psi. */
	double cosi,sini,sin22,cos22;
	double temp;
	double ux = Photon_Ptr->ux;
	double uy = Photon_Ptr->uy;
	double uz = Photon_Ptr->uz;
	double psi;
	double t,I0,I;
	double s11[1000];
	double s12[1000];
	double s33[1000];
	double s43[1000];
	double S[4];
	double S2[4];

	int ithedeg;
	int j,i;

	//cost = SpinTheta(g);
	//sint = sqrt(1.0 - cost*cost);
	/* sqrt() is faster than sin(). */


	// Different layers/object have different scatterers, calculated for skin layers.

	if(Photon_Ptr->inObj==0)
		{
		for(i=0; i<1000;i++)
		  {
		s11[i]=Stokes_Ptr->s11[Photon_Ptr->layer][i];
		s12[i]=Stokes_Ptr->s12[Photon_Ptr->layer][i];
		s33[i]=Stokes_Ptr->s33[Photon_Ptr->layer][i];
		s43[i]=Stokes_Ptr->s43[Photon_Ptr->layer][i];
		}
		}
		else{
			for(i=0; i<1000;i++)
				  {
				s11[i]=Stokes_Ptr->s11[0][i];
				s12[i]=Stokes_Ptr->s12[0][i];
				s33[i]=Stokes_Ptr->s33[0][i];
				s43[i]=Stokes_Ptr->s43[0][i];
				}
		}
	for (j=0;j<4;j++)
	{
		S[j]=Stokes_Ptr->S[j];
	}
	for (j=0;j<4;j++)
	{
	S2[j]=Stokes_Ptr->S2[j];
	}
	//psi = 2.0*PI*RandomNum; /* spin psi 0-2pi. */

	double abc;

	//REJECTION

	do{ t = acos(2*RandomNum1-1);

	  			    psi = RandomNum1*2.0*PI;

	  				I0=s11[0]*S[0]+s12[0]*(S[1]*cos(2*psi)+S[2]*sin(2*psi));

	  	 			ithedeg = floor(t*1000/PI);

	  				I=s11[ithedeg]*S[0]+s12[ithedeg]*(S[1]*cos(2*psi)+S[2]*sin(2*psi));

	  				abc=RandomNum1;

	  			}while(abc*I0>=I);

	cost=cos(t);
	sint=sqrt(1.0-cost*cost);

	//psi = 2.0*PI*RandomNum; /* spin psi 0-2pi. */
	cosp = cos(psi);
	if(psi<PI)
		sinp = sqrt(1.0 - cosp*cosp);
	/* sqrt() is faster than sin(). */
	else
		sinp = - sqrt(1.0 - cosp*cosp);
	//UPDATE Direction Cosines.
	if(fabs(uz) > COSZERO)  { 	/* normal incident. */
		Photon_Ptr->ux = sint*cosp;
		Photon_Ptr->uy = sint*sinp;
		Photon_Ptr->uz = cost*SIGN(uz);
		/* SIGN() is faster than division. */
	}
	else  {		/* regular incident. */
		 temp = sqrt(1.0 - uz*uz);
		Photon_Ptr->ux = sint*(ux*uz*cosp - uy*sinp)
			/temp + ux*cost;
		Photon_Ptr->uy = sint*(uy*uz*cosp + ux*sinp)
			/temp + uy*cost;
		Photon_Ptr->uz = -sint*cosp*temp + uz*cost;
	}

	rotSphi(S, psi, S2);

	S[0]= s11[ithedeg]*S2[0]+s12[ithedeg]*S2[1];

	  			S[1]= s12[ithedeg]*S2[0]+s11[ithedeg]*S2[1];

	  			S[2]= s33[ithedeg]*S2[2]+s43[ithedeg]*S2[3];

	  			S[3]= -s43[ithedeg]*S2[2]+s33[ithedeg]*S2[3];


	  			 temp=(sqrt(1-cost*cost)*sqrt(1-Photon_Ptr->uz*Photon_Ptr->uz));

	  			if ( temp==0){
	  				cosi=0;}
	  			else{

	  				if ((psi>PI) & (psi<2*PI))
	  					cosi=(Photon_Ptr->uz*cost-uz)/temp;
	  				else
	  					cosi=-(Photon_Ptr->uz*cost-uz)/temp;
	  			if (cosi>1) cosi=1;
	  			if (cosi<-1) cosi=-1;
	  			}

	  			sini = sqrt(1-cosi*cosi);

	  			cos22=2*cosi*cosi-1;

	  			sin22=2*sini*cosi;

	  			S2[0]=S[0];

	  			S2[1]=(S[1]*cos22-S[2]*sin22);

	  			S2[2]=(S[1]*sin22+S[2]*cos22);

	  			S2[3]=S[3];


	  			S[1]= S2[1]/S2[0];
	  			S[2]= S2[2]/S2[0];
	  			S[3]= S2[3]/S2[0];
	  			S[0]= 1.0;

	  			for ( j=0;j<4;j++)
	  			{Stokes_Ptr->S[j]=S[j];
	  			}
	  			for (j=0;j<4;j++)
	  			{Stokes_Ptr->S2[j]=S2[j];
	  			}
}

/***********************************************************
*	Move the photon s away in the current layer of medium.
****/
void Hop(PhotonStruct *	Photon_Ptr)
{
	double s = Photon_Ptr->s;

	Photon_Ptr->x += s*Photon_Ptr->ux;
	Photon_Ptr->y += s*Photon_Ptr->uy;
	Photon_Ptr->z += s*Photon_Ptr->uz;
	//if(Photon_Ptr->z<0) Photon_Ptr->z=0;

	Photon_Ptr->pathlength[Photon_Ptr->n_photon-Photon_Ptr->i_photon]+=s;					// Calculates the total pathlength of each photon.
	if(Photon_Ptr->inObj==1)
		Photon_Ptr->pathlengthinobj[Photon_Ptr->n_photon-Photon_Ptr->i_photon]+=s;			// Calculates the total pathlength of each photon in each object.
}

/***********************************************************
*	If uz != 0, return the photon step size in glass,
*	Otherwise, return 0.
*
*	The step size is the distance between the current
*	position and the boundary in the photon direction.
*
*	Make sure uz !=0 before calling this function.
****/
void StepSizeInGlass(PhotonStruct *  Photon_Ptr,
	InputStruct  *  In_Ptr)
{
	double dl_b;	/* step size to boundary. */
	short  layer = Photon_Ptr->layer;
	double uz = Photon_Ptr->uz;

	/* Stepsize to the boundary. */
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z)
		/uz;
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z)
		/uz;
	else
		dl_b = 0.0;

	Photon_Ptr->s = dl_b;
}

/***********************************************************
*	Pick a step size for a photon packet when it is in
*	tissue.
*	If the member sleft is zero, make a new step size
*	with: -log(rnd)/(mua+mus).
*	Otherwise, pick up the leftover in sleft.
*
*	Layer is the index to layer.
*	In_Ptr is the input parameters.
****/
void StepSizeInTissue(PhotonStruct * Photon_Ptr,
	InputStruct  * In_Ptr)
{
	short  layer = Photon_Ptr->layer;
	double mua;
	double mus;
	mua = In_Ptr->layerspecs[layer].mua;
	mus = In_Ptr->layerspecs[layer].mus;

	if(Photon_Ptr->sleft == 0.0) {  /* make a new step. */
		double rnd;

		do rnd = RandomNum1;
		while( rnd <= 0.0 );    /* avoid zero. */
		Photon_Ptr->s = -log(rnd)/(mua+mus);
	}
	else {	/* take the leftover. */
		Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
		Photon_Ptr->sleft = 0.0;
	}

}

/***********************************************************
**	Pick a step size for a photon packet when it is in object.
**	If the member sleft is zero, make a new step size
**	with: -log(rnd)/(mua+mus).
**	Otherwise, pick up the leftover in sleft.
*
**	In_Ptr is the input parameters.
****/
void StepSizeInObject(PhotonStruct * Photon_Ptr,
	InputStruct  * In_Ptr)
{
	double mua;
	double mus;
	mua = In_Ptr->ObjSpecs[0].mua;
	mus = In_Ptr->ObjSpecs[0].mus;

	if(Photon_Ptr->sleft == 0.0) {  /** make a new step. **/
		double rnd;

		do rnd = RandomNum();
		while( rnd <= 0.0 );    /** avoid zero. **/
		Photon_Ptr->s = -log(rnd)/(mua+mus);
	}
	else {	/** take the leftover. **/
		Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
		Photon_Ptr->sleft = 0.0;
	}
}


/***********************************************************
*	Check if the step will hit the boundary.
*	Return 1 if hit boundary.
*	Return 0 otherwise.
*
* 	If the projected step hits the boundary, the members
*	s and sleft of Photon_Ptr are updated.
****/
Boolean HitBoundary(PhotonStruct *  Photon_Ptr,
					InputStruct  *  In_Ptr)
{
	double dl_b;  /* length to boundary. */
	short  layer = Photon_Ptr->layer;
	double uz = Photon_Ptr->uz;
	Boolean hit;


	/* Distance to the boundary. */
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1
				- Photon_Ptr->z)/uz;	/* dl_b>0. */
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0
				- Photon_Ptr->z)/uz;	/* dl_b>0. */

	if(uz != 0.0 && Photon_Ptr->s > dl_b) {
		/* not horizontal & crossing. */
		double mut = In_Ptr->layerspecs[layer].mua
					 + In_Ptr->layerspecs[layer].mus;

		Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
		Photon_Ptr->s     = dl_b;      //Removed temporarily, needs to be there for multi-layer.
		hit = 1;
	}
	else
		hit = 0;

	return(hit);
}

/************************************
** Check if the photon hits the sphere boundary
*******/

Boolean HitSph(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, cr;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double X1, Y1, Z1, X2, Y2, Z2;
	double dB1, dB2;
	double alpha1 = 0, alpha2 = 0.0, beta1 = 0, beta2 = 0.0, gamma1 = 0.0, gamma2 = 0.0;
	double u1, v1, w1, mag1, u2, v2, w2, mag2;
	double mut;
	double err =  1E-6;		/** Very small error **/
	double large = 100000;	/** Large number**/

	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm **/

	A = ux*ux + uy*uy + uz*uz;
	B = 2*(ux*(rx-cx) + uy*(ry-cy) + uz*(rz - cz));
	C = (rx-cx)*(rx-cx) + (ry-cy)*(ry-cy) + (rz - cz)*(rz - cz) - cr*cr;

	del = B*B - 4*A*C;
	if (del < err) { /** tangential case or not hitting case **/
		return(0);
	}
	else { /** If delta is less than 0 then photon path is not intersecting the sphere **/

		srDel = sqrt(del);
		rot1 = (-B + srDel)/(2*A);
		rot2 = (-B - srDel)/(2*A);

		X1 = rx + rot1*ux;
		Y1 = ry + rot1*uy;
		Z1 = rz + rot1*uz;

		X2 = rx + rot2*ux;
		Y2 = ry + rot2*uy;
		Z2 = rz + rot2*uz;

		u1 = X1 - rx;
		v1 = Y1 - ry;
		w1 = Z1 - rz;
		mag1 = sqrt(u1*u1 + v1*v1 + w1*w1);
		if (mag1 != 0)	{ /** If the photon is not at the intersection1, find angle and distance between photon and intersection1 **/
			alpha1 = u1/mag1;
			beta1 = v1/mag1;
			gamma1 = w1/mag1;
			dB1 = (X1-rx)*(X1-rx) + (Y1-ry)*(Y1-ry) + (Z1-rz)*(Z1-rz);
		}
		else	{
			dB1 = large;
		}

		u2 = X2 - rx;
		v2 = Y2 - ry;
		w2 = Z2 - rz;
		mag2 = sqrt(u2*u2 + v2*v2 + w2*w2);
		if (mag2 != 0)	{ /** If the photon is not at the intersection2, find angle and distance between photon and intersection2 **/
			alpha2 = u2/mag2;
			beta2 = v2/mag2;
			gamma2 = w2/mag2;
			dB2 = (X2-rx)*(X2-rx) + (Y2-ry)*(Y2-ry) + (Z2-rz)*(Z2-rz);
		}
		else
		{
			dB2 = large;
		}

		if (fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err
			&& fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** Both the intersection points are in the direction of photon **/
		{
			if (Photon_Ptr->inObj == 0) dl_b = sqrt(dB1<dB2?dB1:dB2); /** For photon outside object **/
			else dl_b = sqrt(dB1>dB2?dB1:dB2); /** For small difference in the boundary and photon location **/
		}
		else if (mag1 !=0 && fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err) /** photon within the sphere **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB1);
		}
		else if (mag2 !=0 && fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** photon within the sphere **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB2);
		}
	}

	if(dl_b == -1)
		return(0);
	else {
		if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
			if (Photon_Ptr->inObj == 1)
				mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
			else
				mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

			Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
			Photon_Ptr->s = dl_b;
			return(1);
		}
		else return(0);
	}
}


/************************************
** Check if the photon hits the cylinder boundary
*******/

Boolean HitCyl(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, cr, cux, cuy, cuz;
	double d, e_x, e_y, e_z, f, g_x, g_y, g_z, h_x, h_y, h_z, delP_x, delP_y, delP_z;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double X1, Y1, Z1, X2, Y2, Z2;
	double dB1, dB2;
	double alpha1 = 0, alpha2 = 0.0, beta1 = 0, beta2 = 0.0, gamma1 = 0.0, gamma2 = 0.0;
	double u1, v1, w1, mag1, u2, v2, w2, mag2;
	double mut;
	double err =  1E-6;		/** Very small error **/
	double large = 1E24;	/** Large number**/

	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	cux = 1;
	cuy = 0;
	cuz = 0;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference  http://www.mrl.nyu.edu/~dzorin/rendering/lectures/lecture3/lecture3-6pp.pdf **/
	d = ux*cux + uy*cuy + uz*cuz;
	e_x = ux-d*cux;
	e_y = uy-d*cuy;
	e_z = uz-d*cuz;

	A = e_x*e_x + e_y*e_y + e_z*e_z ;

	delP_x = rx-cx;
	delP_y = ry-cy;
	delP_z = rz-cz;

	f = delP_x*cux + delP_y*cuy + delP_z*cuz;
	g_x = delP_x - f*cux;
	g_y = delP_y - f*cuy;
	g_z = delP_z - f*cuz;

	B = 2*(e_x*g_x + e_y*g_y + e_z*g_z);

	C = g_x*g_x + g_y*g_y + g_z*g_z - cr*cr;

	del = B*B - 4*A*C;
	if (del < err) { /** tangential case or not hitting case **/
		return(0);
	}
	else { /** If delta is less than 0 then photon path is not intersecting the object **/

		srDel = sqrt(del);
		rot1 = (-B + srDel)/(2*A);
		rot2 = (-B - srDel)/(2*A);

		X1 = rx + rot1*ux;
		Y1 = ry + rot1*uy;
		Z1 = rz + rot1*uz;

		X2 = rx + rot2*ux;
		Y2 = ry + rot2*uy;
		Z2 = rz + rot2*uz;

		u1 = X1 - rx;
		v1 = Y1 - ry;
		w1 = Z1 - rz;
		mag1 = sqrt(u1*u1 + v1*v1 + w1*w1);
		if (mag1 != 0)	{ /** If the photon is not at the intersection1, find angle and distance between photon and intersection1 **/
			alpha1 = u1/mag1;
			beta1 = v1/mag1;
			gamma1 = w1/mag1;
			dB1 = (X1-rx)*(X1-rx) + (Y1-ry)*(Y1-ry) + (Z1-rz)*(Z1-rz);
		}
		else	{
			dB1 = large;
		}

		u2 = X2 - rx;
		v2 = Y2 - ry;
		w2 = Z2 - rz;
		mag2 = sqrt(u2*u2 + v2*v2 + w2*w2);
		if (mag2 != 0)	{ /** If the photon is not at the intersection2, find angle and distance between photon and intersection2 **/
			alpha2 = u2/mag2;
			beta2 = v2/mag2;
			gamma2 = w2/mag2;
			dB2 = (X2-rx)*(X2-rx) + (Y2-ry)*(Y2-ry) + (Z2-rz)*(Z2-rz);
		}
		else
		{
			dB2 = large;
		}

		if (fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err
			&& fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** Both the intersection points are in the direction of photon **/
		{
			if (Photon_Ptr->inObj == 0) dl_b = sqrt(dB1<dB2?dB1:dB2); /** For photon outside object **/
			else dl_b = sqrt(dB1>dB2?dB1:dB2); /** For small difference in the boundary and photon location **/
		}
		else if (mag1 !=0 && fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err) /** photon within the cylinder **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB1);
		}
		else if (mag2 !=0 && fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** photon within the cylinder **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB2);
		}
	}

	if(dl_b == -1)
		return(0);
	else {
		if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
			if (Photon_Ptr->inObj == 1)
				mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
			else
				mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

			Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
			Photon_Ptr->s = dl_b;
			return(1);
		}
		else return(0);
	}
}

/************************************
** Check if the photon hits the ellipsoid boundary
*******/

Boolean HitEll(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	double dl_b = -1;  /** length to boundary. **/
	short  layer = Photon_Ptr->layer;
	double rx, ry, rz;
	double ux, uy, uz;
	double cx, cy, cz, erx, ery, erz;
	double A, B, C;
	double del, srDel, rot1, rot2;
	double X1, Y1, Z1, X2, Y2, Z2;
	double dB1, dB2;
	double alpha1 = 0, alpha2 = 0.0, beta1 = 0, beta2 = 0.0, gamma1 = 0.0, gamma2 = 0.0;
	double u1, v1, w1, mag1, u2, v2, w2, mag2;
	double mut;
	double err =  1E-6;		/** Very small error **/
	double large = 100000;	/** Large number**/

	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	erx = In_Ptr->ObjSpecs[0].crx;
	ery = In_Ptr->ObjSpecs[0].cry;
	erz = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Reference http://www.ogre3d.org/forums/viewtopic.php?f=2&t=26442 **/

	A = (ux*ux)/(erx*erx) + (uy*uy)/(ery*ery) + (uz*uz)/(erz*erz);
	B = ((2*ux*(rx-cx)/(erx*erx)) + (2*uy*(ry-cy)/(ery*ery)) + (2*uz*(rz - cz)/(erz*erz)));
	C = (rx-cx)*(rx-cx)/(erx*erx) + (ry-cy)*(ry-cy)/(ery*ery) + (rz - cz)*(rz - cz)/(erz*erz) - 1;

	del = B*B - 4*A*C;
	if (del < err) { /** tangential case or not hitting case **/
		return(0);
	}
	else { /** If delta is less than 0 then photon path is not intersecting the object **/

		srDel = sqrt(del);
		rot1 = (-B + srDel)/(2*A);
		rot2 = (-B - srDel)/(2*A);

		X1 = rx + rot1*ux;
		Y1 = ry + rot1*uy;
		Z1 = rz + rot1*uz;

		X2 = rx + rot2*ux;
		Y2 = ry + rot2*uy;
		Z2 = rz + rot2*uz;

		u1 = X1 - rx;
		v1 = Y1 - ry;
		w1 = Z1 - rz;
		mag1 = sqrt(u1*u1 + v1*v1 + w1*w1);
		if (mag1 != 0)	{ /** If the photon is not at the intersection1, find angle and distance between photon and intersection1 **/
			alpha1 = u1/mag1;
			beta1 = v1/mag1;
			gamma1 = w1/mag1;
			dB1 = (X1-rx)*(X1-rx) + (Y1-ry)*(Y1-ry) + (Z1-rz)*(Z1-rz);
		}
		else {
			dB1 = large;
		}

		u2 = X2 - rx;
		v2 = Y2 - ry;
		w2 = Z2 - rz;
		mag2 = sqrt(u2*u2 + v2*v2 + w2*w2);
		if (mag2 != 0)	{ /** If the photon is not at the intersection2, find angle and distance between photon and intersection2 **/
			alpha2 = u2/mag2;
			beta2 = v2/mag2;
			gamma2 = w2/mag2;
			dB2 = (X2-rx)*(X2-rx) + (Y2-ry)*(Y2-ry) + (Z2-rz)*(Z2-rz);
		}
		else
		{
			dB2 = large;
		}

		if (fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err
			&& fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** Both the intersection points are in the direction of photon **/
		{
			if (Photon_Ptr->inObj == 0) dl_b = sqrt(dB1<dB2?dB1:dB2); /** For photon outside ellipsoid **/
			else dl_b = sqrt(dB1>dB2?dB1:dB2); /** For small difference in the boundary and photon location **/
		}
		else if (mag1 !=0 && fabs(alpha1-ux) < err && fabs(beta1-uy) < err && fabs(gamma1-uz) < err) /** photon within the ellipsoid **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB1);
		}
		else if (mag2 !=0 && fabs(alpha2-ux) < err && fabs(beta2-uy) < err && fabs(gamma2-uz) < err) /** photon within the ellipsoid **/
		{
			if (Photon_Ptr->inObj == 0) return (0);
			else dl_b = sqrt(dB2);
		}
	}

	if(dl_b == -1)
		return(0);
	else {
		if(Photon_Ptr->s >= dl_b) { /** Check with stepsize if the photon reaches the boundary **/
			if (Photon_Ptr->inObj == 1)
				mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
			else
				mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;

			Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
			Photon_Ptr->s = dl_b;
			return(1);
		}
		else return(0);
	}
}

/***********************************************************
**	Check if the given distance moved by a photon
**	lies on the surface of box.
**	Return 1 if hits surface.
**	Return 0 otherwise.
*****/

Boolean checkOnPlane(double bXMin, double bYMin, double bZMin, double bXMax, double bYMax, double bZMax,
					 double x, double y, double z, double ux, double uy, double uz, double dist)
{
	double newX, newY, newZ;
	newX = (x + ux*dist);	/*Find the point after the photon travels the distance */
	newY = (y + uy*dist);
	newZ = (z + uz*dist);
	if((bXMin <= newX && newX <= bXMax) && /* Check if the point is on box surface */
	   (bYMin <= newY && newY <= bYMax) &&
	   (bZMin <= fabs(newZ) && newZ <= bZMax)) {
		   return(1);
	}
	else return(0);
}
/***********************************************************
**	Check if the step will hit the surface of box.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
**	The minimum distance to boundary is saved.
****/

Boolean hitPla(int cubNum, InputStruct  *  In_Ptr,
			   PhotonStruct *  Photon_Ptr, double * minDistBox)
{
	double cXMin, cYMin, cZMin, cXMax, cYMax, cZMax;
	double x, y, z, ux, uy, uz;
	double distX1, distX2, distY1, distY2, distZ1, distZ2;
	double minDistX = 9999, minDistY = 9999, minDistZ = 9999;
	float newX, newY, newZ;
	double minDist, tempMin;
	double ERR = 1E-10, HIGHVAL = 99999;

	/** Get the box dimensions **/
	cXMin = In_Ptr->ObjSpecs[cubNum].XMin;
	cYMin = In_Ptr->ObjSpecs[cubNum].YMin;
	cZMin = In_Ptr->ObjSpecs[cubNum].ZMin;
	cXMax = In_Ptr->ObjSpecs[cubNum].XMax;
	cYMax = In_Ptr->ObjSpecs[cubNum].YMax;
	cZMax = In_Ptr->ObjSpecs[cubNum].ZMax;

	x = Photon_Ptr->x;
	y = Photon_Ptr->y;
	z = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/**	Find the distance of photon with respect to 6 directions. Discard the negative distance	**/
	if(ux != 0) {
		distX1 = (cXMin-x)/ux;
		if(distX1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distX1)){
			distX1 = HIGHVAL;		/** Setting a high value because point doesn't lie on surface with the distance **/
		}
		distX2 = (cXMax-x)/ux;
		if(distX2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distX2)){
			distX2=HIGHVAL;
		}
		minDistX = distX1<distX2?distX1:distX2;
	}

	if(uy != 0) {
		distY1 = (cYMin-y)/uy;
		if(distY1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distY1)){
			distY1 = HIGHVAL;
		}
		distY2 = (cYMax-y)/uy;
		if(distY2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distY2)){
			distY2 = HIGHVAL;
		}
		minDistY = distY1<distY2?distY1:distY2;
	}

	if (uz !=0) {
		distZ1 = (cZMin-z)/uz;
		if(distZ1<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distZ1)){
			distZ1 = HIGHVAL;
		}
		distZ2 = (cZMax-z)/uz;
		if(distZ2<ERR || !checkOnPlane(cXMin, cYMin, cZMin, cXMax, cYMax, cZMax, x, y, z, ux, uy, uz, distZ2)){
			distZ2 = HIGHVAL;
		}
		minDistZ = distZ1<distZ2?distZ1:distZ2;
	}

	if (minDistX == HIGHVAL && minDistY == HIGHVAL && minDistZ == HIGHVAL) {
		return 0; /** No hit when there is no valid distance. **/
	}
	minDist=(minDistX<minDistY?minDistX:minDistY)<minDistZ?(minDistX<minDistY?minDistX:minDistY):minDistZ;
	//minDist=min(min(minDistX, minDistY), minDistZ); /** Minimum of the distance in x, y, and z directions **/

	if (minDist <= Photon_Ptr->s) {
		*minDistBox = minDist;
		return(1);		/** Hit if minimum distance is less than step-size **/
	}
	return(0);
}

/***********************************************************
**	Check if the step will hit the surface of cuboid.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
*
**	If the projected step hits the boundary, the members
**	s and sleft of Photon_Ptr are updated.
****/
Boolean HitCub(PhotonStruct *  Photon_Ptr,
			   InputStruct  *  In_Ptr)
{
	Boolean hitPlane = 0;
	double hitDist;
	double mut;

	hitPlane = hitPla(0, In_Ptr, Photon_Ptr, &hitDist);		/** Check if photon hits cuboid **/

	/**	Update the photon pointer when there is a hit. **/
	if(hitPlane == 1) {
		if(hitDist <= Photon_Ptr->s) {
			if (Photon_Ptr->inObj == 1) {
				mut = In_Ptr->ObjSpecs[0].mua + In_Ptr->ObjSpecs[0].mus;
			}
			else {
				mut = In_Ptr->layerspecs[Photon_Ptr->layer].mua + In_Ptr->layerspecs[Photon_Ptr->layer].mus;
			}

			Photon_Ptr->sleft = (Photon_Ptr->s - hitDist)*mut;
			Photon_Ptr->s = hitDist;
			return(1);
		}
	}
	else{
		return(0);
	}
}

/***********************************************************
*	Drop photon weight inside the tissue (not glass).
*
*  The photon is assumed not dead.
*
*	The weight drop is dw = w*mua/(mua+mus).
*
*	The dropped weight is assigned to the absorption array
*	elements.
****/
void Drop(InputStruct  *	In_Ptr,
		  PhotonStruct *	Photon_Ptr,
		  Stokes_struct * Stokes_Ptr,
		  OutStruct * Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double izd, ird;	/* LW 5/20/98. To avoid out of short range.*/
	short  ix, iy, iz, ir;	/* index to z & r. */
	short  layer = Photon_Ptr->layer;
	double mua, mus;
	short i;

	/* compute array indices. */
	ix=iy=199;  /*To avoid junk value. All negative indices are also stored in the last bin.*/

	if(x>=-1) ix=(short)((x+1)/In_Ptr->dx);
	if(ix>199) ix=199;

	if(y>=-1) iy=(short)((y+1)/In_Ptr->dy);
	if(iy>199) iy=199;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd>In_Ptr->nz-1) iz=In_Ptr->nz-1;
	else iz = izd;

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	/* update photon weight. */
	mua = In_Ptr->layerspecs[layer].mua;
	mus = In_Ptr->layerspecs[layer].mus;
	dwa = Photon_Ptr->w * mua/(mua+mus);
	Photon_Ptr->w -= dwa;


	/* assign dwa to the absorption array element. */
	//Out_Ptr->A_rz[ir][iz] += dwa;

//	if(ix==100&&iy==100)
	//Out_Ptr->P_z[iz] +=Photon_Ptr->w;


	if(iy==100&&ix>=0)
		Out_Ptr->A_xz[iz][ix]+=dwa;

}

/***********************************************************
**	Drop photon weight inside the object.
*
**  The photon is assumed not dead.
*
**	The weight drop is dw = w*mua/(mua+mus).
*
**	The dropped weight is assigned to the absorption variable
**	elements in rz coordinates and object absorbance variable.
****/
void DropObj(InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 Stokes_struct * Stokes_Ptr,
			 OutStruct *	 Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double izd, ird;	/* LW 5/20/98. To avoid out of short range.*/
	short  ix ,iy ,iz, ir;	/* index to z & r. */
	short  layer = Photon_Ptr->layer;
	double mua, mus;
	short i;

	/* compute array indices. */

	ix=iy=199;  /*To avoid junk value. All negative indices are also stored in the last bin.*/

	if(x>=-1) ix=(short)((x+1)/In_Ptr->dx);
		if(ix>199) ix=199;

	if(y>=-1) iy=(short)((y+1)/In_Ptr->dy);
		if(iy>199) iy=199;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz - 1) iz = In_Ptr->nz - 1;
	else iz = izd;

	ird = sqrt(x*x + y*y)/In_Ptr->dr;
	if(ird > In_Ptr->nr - 1) ir = In_Ptr->nr - 1;
	else ir = ird;

	/** update photon weight. **/
	mua = In_Ptr->ObjSpecs[0].mua;
	mus = In_Ptr->ObjSpecs[0].mus;
	dwa = Photon_Ptr->w * mua/(mua+mus);
	Photon_Ptr->w -= dwa;

//	Out_Ptr->A_rz[ir][iz] += dwa;	/* assign dwa to the absorption array element. */

	if(iy==100)
			Out_Ptr->A_xz[iz][ix]+=dwa;

	//if(ix==100&&iy==100)
	//Out_Ptr->P_z[iz] += Photon_Ptr->w;

	//Out_Ptr->A_obj += dwa;	/** assign dwa to the absorption variable within object. **/



}

/***********************************************************
*	The photon weight is small, and the photon packet tries
*	to survive a roulette.
****/
void Roulette(PhotonStruct * Photon_Ptr,InputStruct * In_Ptr)
{
	double rnd=RandomNum1;
	if( Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead)
	{

	if(Photon_Ptr->w == 0.0)
		Photon_Ptr->dead = 1;
	else if(rnd < CHANCE) /* survived the roulette.*/
		Photon_Ptr->w /= CHANCE;
	else
		Photon_Ptr->dead = 1;
	}
}


/***********************************************************
*	Record the photon weight exiting the first layer(uz<0),
*	no matter whether the layer is glass or not, to the
*	reflection array.
*
*	Update the photon weight as well.
****/
void RecordR(double			Refl,	/* reflectance. */
			 InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 Stokes_struct * Stokes_Ptr,
			 OutStruct	  *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy, ir, ia;	/* index to r & angle. */
	double ird, iad;	/* LW 5/20/98. To avoid out of short range.*/
	double phi;
	//phi=atan2(Photon_Ptr->uy,Photon_Ptr->ux);
	//rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);

	ix=iy=199;  /*To avoid junk value. All negative indices are also stored in the last bin.*/

	if(x>=-0.1) ix=(short)((x+0.1)/(In_Ptr->dx/10));
		if(ix>199) ix=199;

	if(y>=-0.1) iy=(short)((y+0.1)/(In_Ptr->dy/10));
		if(iy>199) iy=199;

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	iad = acos(-Photon_Ptr->uz)/In_Ptr->da;
	if(iad>In_Ptr->na-1) ia=In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the reflection array element. */
	//Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);

	if(ix>=0&&iy>=0)

	{
	Out_Ptr->Rd_xy[iy][ix] += Photon_Ptr->w*(1.0-Refl);

	Out_Ptr->IR[iy][ix] += Stokes_Ptr->S[0];

	Out_Ptr->QR[iy][ix] += Stokes_Ptr->S[1];

	Out_Ptr->UR[iy][ix] += Stokes_Ptr->S[2];

	Out_Ptr->VR[iy][ix] += Stokes_Ptr->S[3];

	Out_Ptr->IR_1+=Stokes_Ptr->S[0];

	Out_Ptr->QR_1+=Stokes_Ptr->S[1];

	Out_Ptr->UR_1+=Stokes_Ptr->S[2];

	Out_Ptr->VR_1+=Stokes_Ptr->S[3];

	}

	Photon_Ptr->w *= Refl;  //removed temporarily, needed for mcmleo
}


/***********************************************************
*	Record the photon weight exiting the last layer(uz>0),
*	no matter whether the layer is glass or not, to the
*	transmittance array.
*
*	Update the photon weight as well.
****/
void RecordT(double 		    Refl,
			 InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 Stokes_struct * Stokes_Ptr,
			 OutStruct    *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy, ir, ia;	/* index to r & angle. */
	double ird, iad;	/* LW 5/20/98. To avoid out of short range.*/
	double phi;
	//phi=-atan2(Photon_Ptr->uy,Photon_Ptr->ux);
    //rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);

	ix=iy=199;  /*To avoid junk value. All negative indices are also stored in the last bin.*/

	if(x>=-1) ix=(short)((x+1)/In_Ptr->dx);
			if(ix>199) ix=199;

	if(y>=-1) iy=(short)((y+1)/In_Ptr->dy);
			if(iy>199) iy=199;

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird>In_Ptr->nr-1) ir=In_Ptr->nr-1;
	else ir = ird;

	iad = acos(Photon_Ptr->uz)/In_Ptr->da; /* LW 1/12/2000. Removed -. */
	if(iad>In_Ptr->na-1) ia=In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the transmittance array element. */
	//Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);


	if(ix>=0&&iy>=0)
	{
	Out_Ptr->Tt_xy[iy][ix] += Photon_Ptr->w*(1.0-Refl);

	Out_Ptr->ITT[iy][ix] += Stokes_Ptr->S[0]*Photon_Ptr->w;

	Out_Ptr->QTT[iy][ix] += Stokes_Ptr->S[1]*Photon_Ptr->w;

	Out_Ptr->UTT[iy][ix] += Stokes_Ptr->S[2]*Photon_Ptr->w;

	Out_Ptr->VTT[iy][ix] += Stokes_Ptr->S[3]*Photon_Ptr->w;
	}
	Photon_Ptr->w *= Refl;//removed temporarily, needed for mcmleo

}


/***********************************************************
*	Decide whether the photon will be transmitted or
*	reflected on the upper boundary (uz<0) of the current
*	layer.
*
*	If "layer" is the first layer, the photon packet will
*	be partially transmitted and partially reflected if
*	PARTIALREFLECTION is set to 1,
*	or the photon packet will be either transmitted or
*	reflected determined statistically if PARTIALREFLECTION
*	is set to 0.
*
*	Record the transmitted photon weight as reflection.
*
*	If the "layer" is not the first layer and the photon
*	packet is transmitted, move the photon to "layer-1".
*
*	Update the photon parmameters.
****/
void CrossUpOrNot(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr,
				  Stokes_struct * Stokes_Ptr,
				  OutStruct	   *Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. always */
	double phi; /*For Rotsphi*/
	short j;
	/* positive. */
	double r=0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer-1].n;
	phi=atan2(Photon_Ptr->uy,Photon_Ptr->ux);

	/* Get r. */
	if( - uz <= In_Ptr->layerspecs[layer].cos_crit0)
		r=1.0;		      /* total internal reflection. */
	else r = RFresnel(Stokes_Ptr,ni, nt, -uz, &uz1);

#if PARTIALREFLECTION
	if(layer == 1 && r<1.0) {	/* partially transmitted. */
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
		StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
		Photon_Ptr->uz = -uz1;	/* transmitted photon. */
		RecordR(r, In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
		Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
								Photon_Ptr,Stokes_Ptr);
		Photon_Ptr->uz = -uz;	/* reflected photon. */
	}
	else if(RandomNum() > r) {/* transmitted to layer-1. */
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
		StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
		Photon_Ptr->layer--;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = -uz1;
		phi=atan2(Photon_Ptr->uy,Photon_Ptr->ux);
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
									{
										Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
									}
	}
	else			      		/* reflected. */
		{
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				StokesReflect(Stokes_Ptr,acos(fabs(uz)),&uz1);
				Photon_Ptr->uz = -uz;
				phi=-phi;
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
									{
										Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
									}
		}
#else
	if(RandomNum() > r) {		/* transmitted to layer-1. */       //very,very important, required for mcmleo
		if(layer==1)  {
			rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
			StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
			Photon_Ptr->uz = -uz1;
			RecordR(0.0, In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
			Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
									Photon_Ptr,Stokes_Ptr);
			Photon_Ptr->dead = 1;
		}
		else {
			rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
			StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
			Photon_Ptr->layer--;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = -uz1;
			phi=atan2(Photon_Ptr->uy,Photon_Ptr->ux);
					rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
					for (j=0;j<4;j++)
										{
											Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
										}
		}
	}																//very,very important, required for mcmleo
	else 						/* reflected. */					//very,very important, required for mcmleo
	{
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				StokesReflect(Stokes_Ptr,acos(fabs(uz)),&uz1);
				Photon_Ptr->uz = -uz;
				phi=-phi;
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
									{
										Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
									}									//very,very important, required for mcmleo
	}
#endif
}



/***********************************************************
*	Decide whether the photon will be transmitted  or be
*	reflected on the bottom boundary (uz>0) of the current
*	layer.
*
*	If the photon is transmitted, move the photon to
*	"layer+1". If "layer" is the last layer, record the
*	transmitted weight as transmittance. See comments for
*	CrossUpOrNot.
*
*	Update the photon parmameters.
****/
void CrossDnOrNot(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr,
				  Stokes_struct * Stokes_Ptr,
				  OutStruct	   *Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. */
	double r=0.0;	/* reflectance */
	double phi;
	short j;
	short  layer = Photon_Ptr->layer;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer+1].n;
	phi=-atan2(Photon_Ptr->uy,Photon_Ptr->ux);
	/* Get r. */
	if( uz <= In_Ptr->layerspecs[layer].cos_crit1)
		r=1.0;		/* total internal reflection. */
	else r = RFresnel(Stokes_Ptr,ni, nt, uz, &uz1);

#if PARTIALREFLECTION
	if(layer == In_Ptr->num_layers && r<1.0) {
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
		StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
		Photon_Ptr->uz = uz1;
		RecordT(r, In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
		Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
								Photon_Ptr,Stokes_Ptr);
		Photon_Ptr->uz = -uz;
	}
	else if(RandomNum() > r) {/* transmitted to layer+1. */
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
		StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
		Photon_Ptr->layer++;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = uz1;
		phi=-atan2(Photon_Ptr->uy,Photon_Ptr->ux);
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
							{
								Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
							}
	}
	else 						/* reflected. */
	{
		rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				StokesReflect(Stokes_Ptr,acos(fabs(uz)),&uz1);
				Photon_Ptr->uz = -uz;
				phi=-phi;
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
									{
										Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
									}
	}
#else
	if(RandomNum() > r) {		/* transmitted to layer+1. */					//very,very important, required for mcmleo
		if(layer == In_Ptr->num_layers) {
			rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
			StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
			Photon_Ptr->uz = uz1;
			RecordT(0.0, In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
			Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
									Photon_Ptr,Stokes_Ptr);
			Photon_Ptr->dead = 1;
			phi=-atan2(Photon_Ptr->uy,Photon_Ptr->ux);
					rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
					for (j=0;j<4;j++)
								{
									Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
								}
		}
		else {
			rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
			StokesTransmit(Stokes_Ptr,acos(fabs(uz)),&uz1,ni,nt);
			Photon_Ptr->layer++;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = uz1;
		}
	}																			//very,very important, required for mcmleo
	else 						/* reflected. */								//very,very important, required for mcmleo
		{rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				StokesReflect(Stokes_Ptr,acos(fabs(uz)),&uz1);
				Photon_Ptr->uz = -uz;
				phi=-phi;
				rotSphi(Stokes_Ptr->S,phi,Stokes_Ptr->S2);
				for (j=0;j<4;j++)
									{
										Stokes_Ptr->S[j]=Stokes_Ptr->S2[j];
									}													//very,very important, required for mcmleo
		}
#endif
}

/***********************************************************
* If uz>0 check for boundary below and if uz<0 check for boundary above
****/
void CrossOrNot(InputStruct  *	In_Ptr,
			    PhotonStruct *	Photon_Ptr,
			    Stokes_struct * Stokes_Ptr,
				OutStruct    *	Out_Ptr)
{
	if(Photon_Ptr->uz < 0.0)
		CrossUpOrNot(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
	else
		CrossDnOrNot(In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
}

/***********************************************************
**	Decide whether the photon will be transmitted or
**	reflected from sphere to layer or vise-versa
*
**	The photon packet will be either transmitted or
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotSph(InputStruct  *	In_Ptr,
				   PhotonStruct *	Photon_Ptr,
				   Stokes_struct * Stokes_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni;
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, cr;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;


	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Change the xyz coordinate system to local coordinate system whose z axis align along the surface normal at the point of intersection **/
	if (Photon_Ptr->inObj == 0){
		u1 = cx-rx;
		v1 = cy-ry;
		w1 = cz-rz;
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = rx-cx;
		v1 = ry-cy;
		w1 = rz-cz;
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;

	if(fabs(gamma1) < COSZERO){
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit)
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(Stokes_Ptr,ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}

/***********************************************************
**	Decide whether the photon will be transmitted or
**	reflected from cylinder to layer or vise-versa
*
**	The photon packet will be either transmitted or
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotCyl(InputStruct  *	In_Ptr,
				   PhotonStruct *	Photon_Ptr,
				   Stokes_struct * Stokes_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni;
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, cr;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;
	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	cr = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;
	if (Photon_Ptr->inObj == 0){
		u1 = rx-rx;
		v1 = cy-ry;
		w1 = cz-rz;
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = rx-rx;
		v1 = ry-cy;
		w1 = rz-cz;
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;

	if(fabs(gamma1) < COSZERO){
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit)
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(Stokes_Ptr,ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}

/***********************************************************
**	Decide whether the photon will be transmitted or
**	reflected from ellipsoid to layer or vise-versa
*
**	The photon packet will be either transmitted or
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotEll(InputStruct  *	In_Ptr,
				   PhotonStruct *	Photon_Ptr,
				   Stokes_struct * Stokes_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni;
	double nt;
	double cos_crit;

	double rx, ry, rz;
	double ux, uy, uz, ux1, uy1, uz1, uz2;
	double cx, cy, cz, erx, ery, erz, d;
	double alpha1 = 0.0, beta1 = 0.0, gamma1 = 0.0;
	double u1, v1, w1, mag1;


	cx = In_Ptr->ObjSpecs[0].cx;
	cy = In_Ptr->ObjSpecs[0].cy;
	cz = In_Ptr->ObjSpecs[0].cz;
	erx = In_Ptr->ObjSpecs[0].crx;
	ery = In_Ptr->ObjSpecs[0].cry;
	erz = In_Ptr->ObjSpecs[0].crz;
	rx = Photon_Ptr->x;
	ry = Photon_Ptr->y;
	rz = Photon_Ptr->z;
	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/** Change the xyz coordinate system to local coordinate system whose z axis align along the surface normal at the point of intersection **/
	if (Photon_Ptr->inObj == 0){
		u1 = ((rx-cx)/(erx*erx));
		v1 = ((ry-cy)/(ery*ery));
		w1 = ((rz-cz)/(erz*erz));
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit0;
	}
	else {
		u1 = ((rx-cx)/(erx*erx));
		v1 = ((ry-cy)/(ery*ery));
		w1 = ((rz-cz)/(erz*erz));
		ni = In_Ptr->ObjSpecs[0].n;
		nt = In_Ptr->layerspecs[layer].n;
		cos_crit = In_Ptr->ObjSpecs[0].cos_crit1;
	}

	mag1   = sqrt(u1*u1 + v1*v1 + w1*w1);
	alpha1 = u1/mag1;
	beta1  = v1/mag1;
	gamma1 = w1/mag1;

	if(fabs(gamma1) < COSZERO){
		ux1 = ux*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uz*(-sqrt(1-gamma1*gamma1));
		uy1 = ux*(-beta1/sqrt(1-gamma1*gamma1)) + uy*(alpha1/sqrt(1-gamma1*gamma1));
		uz1 = ux*alpha1 + uy*beta1 + uz*gamma1;
	}
	else {
		ux1 = ux;
		uy1 = uy;
		uz1 = SIGN(gamma1)*uz;
	}

	/** Get r. **/
	/** Check if the angle of incidence is greater or lesser than critical angle and find the transmission angle **/
	if( fabs(uz1) <= cos_crit)
		r=1.0;		      /** total internal reflection. **/
	else r = RFresnel(Stokes_Ptr,ni, nt, fabs(uz1), &uz2);

	if(RandomNum() > r) {		/** transmitted. **/
		ux1 *= ni/nt;
		uy1 *= ni/nt;
		uz1  = SIGN(uz1)*uz2;
		Photon_Ptr->inObj = (!Photon_Ptr->inObj);
	}
	else {					/** reflected. **/
		uz1 = -uz1;
	}
	/** Convert back to xyz coordinate system **/
	if(fabs(gamma1) < COSZERO){
		Photon_Ptr->ux = ux1*(alpha1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(-beta1/sqrt(1-gamma1*gamma1)) + uz1*alpha1;
		Photon_Ptr->uy = ux1*(beta1*gamma1/sqrt(1-gamma1*gamma1)) + uy1*(alpha1/sqrt(1-gamma1*gamma1)) + uz1*beta1;
		Photon_Ptr->uz = ux1*(-sqrt(1-gamma1*gamma1)) + uz1*gamma1;
	}
	else {
		Photon_Ptr->ux = ux1;
		Photon_Ptr->uy = uy1;
		Photon_Ptr->uz = SIGN(gamma1)*uz1;
	}
}


/***********************************************************
**	Decide whether the photon will be transmitted or
**	reflected from cuboid to layer or vise-versa
*
**	The photon packet will be either transmitted or
**	reflected determined statistically.
*
*	Update the photon parmameters.
****/
void CrossOrNotCub(InputStruct  *	In_Ptr,
				   PhotonStruct *	Photon_Ptr,
				   Stokes_struct * Stokes_Ptr,
				   OutStruct    *	Out_Ptr)
{
	double ux = Photon_Ptr->ux; /* x directional cosine. */
	double uy = Photon_Ptr->uy; /* y directional cosine. */
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double plaAng, critAng, uz1;	/* cosines of transmission alpha. */
	int plane;
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->layer;
	double ni;
	double nt;

	if (Photon_Ptr->inObj == 0) {
		ni = In_Ptr->layerspecs[layer].n;
		nt = In_Ptr->ObjSpecs[0].n;
	}
	else {
		nt = In_Ptr->layerspecs[layer].n;
		ni = In_Ptr->ObjSpecs[0].n;
	}
	if(Photon_Ptr->x == In_Ptr->ObjSpecs[0].XMin) {
		plaAng = Photon_Ptr->ux;
		plane = 1;
	}
	else if(Photon_Ptr->x == In_Ptr->ObjSpecs[0].XMax) {
		plaAng = Photon_Ptr->ux;
		plane = 2;
	}
	else if(Photon_Ptr->y == In_Ptr->ObjSpecs[0].YMin) {
		plaAng = Photon_Ptr->uy;
		plane = 3;
	}
	else if(Photon_Ptr->y == In_Ptr->ObjSpecs[0].YMax) {
		plaAng = Photon_Ptr->uy;
		plane = 4;
	}
	else if(Photon_Ptr->z == In_Ptr->ObjSpecs[0].ZMin) {
		plaAng = Photon_Ptr->uz;
		plane = 5;
	}
	else {
		plaAng = Photon_Ptr->uz;
		plane = 6;
	}

	if (Photon_Ptr->inObj == 1) {
		critAng = In_Ptr->ObjSpecs[0].cos_crit1;
	}
	else {
		critAng = In_Ptr->ObjSpecs[0].cos_crit0;
	}

	/* Get r. */

	if( fabs(plaAng) <= critAng)
		r = 1.0;		/* total internal reflection. */
	else r = RFresnel(Stokes_Ptr,ni, nt, fabs(plaAng), &uz1);

	if(RandomNum() > r) {
		if (Photon_Ptr->inObj == 1)	 {	/*If photon is in material change to container*/
			Photon_Ptr->inObj = 0;
		}
		else {
			Photon_Ptr->inObj = 1;
		}
		switch (plane)
		{
		case 1:
			Photon_Ptr->ux = SIGN(ux)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 2:
			Photon_Ptr->ux = SIGN(ux)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 3:
			Photon_Ptr->uy = SIGN(uy)*uz1;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 4:
			Photon_Ptr->uy = SIGN(uy)*uz1;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uz *= ni/nt;
			break;

		case 5:
			Photon_Ptr->uz = SIGN(uz)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->ux *= ni/nt;
			break;

		case 6:
			Photon_Ptr->uz = SIGN(uz)*uz1;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->ux *= ni/nt;
			break;

		default:
			break;
		}
	}
	else {
		switch (plane)
		{
		case 1:
			Photon_Ptr->ux = -ux;
			break;

		case 2:
			Photon_Ptr->ux = -ux;
			break;

		case 3:
			Photon_Ptr->uy = -uy;
			break;

		case 4:
			Photon_Ptr->uy = -uy;
			break;

		case 5:
			Photon_Ptr->uz = -uz;
			break;

		case 6:
			Photon_Ptr->uz = -uz;
			break;

		default:
			break;
		}
	}
}


/******************
** Check which object the CrossOrNot should check for
*******************/
void CrossOrNotObj(InputStruct  *  In_Ptr,
				   PhotonStruct *  Photon_Ptr,
				   Stokes_struct * Stokes_Ptr,
				   OutStruct    *  Out_Ptr)
{
	switch(In_Ptr->objCode) {
	case 0:
		CrossOrNot(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		break;
	case 1:
		CrossOrNotSph(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		break;
	case 2:
		CrossOrNotCyl(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		break;
	case 3:
		CrossOrNotEll(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		break;
	case 4:
		CrossOrNotCub(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		break;
	}
}

/***********************************************************
*	Move the photon packet in glass layer.
*	Horizontal photons are killed because they will
*	never interact with tissue again.
****/
void HopInGlass(InputStruct  * In_Ptr,
				PhotonStruct * Photon_Ptr,
				Stokes_struct * Stokes_Ptr,
				OutStruct    * Out_Ptr)
{
	if(Photon_Ptr->uz == 0.0) {
		/* horizontal photon in glass is killed. */
		Photon_Ptr->dead = 1;
	}
	else {
		StepSizeInGlass(Photon_Ptr, In_Ptr);
		Hop(Photon_Ptr);
		CrossOrNot(In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
	}
}

/***********************************************************
*	Set a step size, move the photon, drop some weight,
*	choose a new photon direction for propagation.
*
*	When a step size is long enough for the photon to
*	hit an interface, this step is divided into two steps.
*	First, move the photon to the boundary free of
*	absorption or scattering, then decide whether the
*	photon is reflected or transmitted.
*	Then move the photon in the current or transmission
*	medium with the unfinished stepsize to interaction
*	site.  If the unfinished stepsize is still too long,
*	repeat the above process.
****/
void HopDropSpinInTissue(InputStruct  *  In_Ptr,
						 PhotonStruct *  Photon_Ptr,
						 Stokes_struct * Stokes_Ptr,
						 OutStruct    *  Out_Ptr)
{
	StepSizeInTissue(Photon_Ptr, In_Ptr);

	if(HitBoundary(Photon_Ptr, In_Ptr)) {
		Hop(Photon_Ptr);	/* move to boundary plane. */
		//Drop(In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);       //to be removed later for mcmleo
		CrossOrNot(In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
	}
	else {
		Hop(Photon_Ptr);
		Drop(In_Ptr, Photon_Ptr, Stokes_Ptr, Out_Ptr);
		Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
			Photon_Ptr, Stokes_Ptr);
	}
}


/******************
** Check which object the Hit has to be should checked for
*******************/
Boolean HitObj(PhotonStruct * Photon_Ptr,
			   InputStruct  * In_Ptr)
{

	Boolean hit;
	switch (In_Ptr->objCode) {
	case 0:
		hit=HitBoundary(Photon_Ptr, In_Ptr);
		break;
	case 1:
		hit=HitSph(Photon_Ptr, In_Ptr);
		break;
	case 2:
		hit=HitCyl(Photon_Ptr, In_Ptr);
		break;
	case 3:
		hit=HitEll(Photon_Ptr, In_Ptr);
		break;
	case 4:
		hit=HitCub(Photon_Ptr, In_Ptr);
		break;
	}
	return(hit);
}

/***********************************************************
*	Set a step size, move the photon, drop some weight,
*	choose a new photon direction for propagation.
*
*	When a step size is long enough for the photon to
**	hit an interface or object, this step is divided into two steps.
*	First, move the photon to the boundary free of
*	absorption or scattering, then decide whether the
*	photon is reflected or transmitted.
*	Then move the photon in the current or transmission
*	medium with the unfinished stepsize to interaction
*	site.  If the unfinished stepsize is still too long,
*	repeat the above process.
****/
void HopDropSpinInTissueObj(InputStruct  *  In_Ptr,
							PhotonStruct *  Photon_Ptr,
							Stokes_struct * Stokes_Ptr,
							OutStruct    *  Out_Ptr)
{
	if (Photon_Ptr->inObj == 1)
		StepSizeInObject(Photon_Ptr, In_Ptr); /* Stepsize by mut of object */
	else
		StepSizeInTissue(Photon_Ptr, In_Ptr); /* Stepsize by mut of Photon->layer  */

	if (Photon_Ptr->layer == In_Ptr->objLayer){ /** Check if the photon is in the layer where the object is **/
		if (HitObj(Photon_Ptr, In_Ptr)){ /** Check if photon hits object boundary **/
			Hop(Photon_Ptr);	/* move to object surface. */
			CrossOrNotObj(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		}
		else { /** Photon is in the layer of object and doesn't hit object **/
			if(Photon_Ptr->inObj == 0) { /** Check if photon hits layer boundary if it is outside object **/
				if(HitBoundary(Photon_Ptr, In_Ptr)) { /* If photon hits layer boundary */
					Hop(Photon_Ptr);	/* move to boundary plane. */
					CrossOrNot(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
				}
				else { /* Hop,drop,spin in layer */
					Hop(Photon_Ptr);
					Drop(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
					Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
						Photon_Ptr,Stokes_Ptr);
				}
			}
			else { /** Photon in object. Hop, drop, spin in sphere.**/
				Hop(Photon_Ptr);
				DropObj(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
				Spin(In_Ptr->ObjSpecs[0].g, Photon_Ptr,Stokes_Ptr);
			}
		}
	}
	else { /** Photon not in the layer of object. Conventional boundary check and hop, drop, spin **/
		if(HitBoundary(Photon_Ptr, In_Ptr)) {
			Hop(Photon_Ptr);	/* move to boundary plane. */
			CrossOrNot(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		}
		else {
			Hop(Photon_Ptr);
			Drop(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
			Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g,
				Photon_Ptr,Stokes_Ptr);
		}
	}
}


/***********************************************************
****/
void HopDropSpin(InputStruct  *  In_Ptr,
				 PhotonStruct *  Photon_Ptr,
				 Stokes_struct * Stokes_Ptr,
				 OutStruct    *  Out_Ptr)
{
	short layer = Photon_Ptr->layer;

	if((In_Ptr->layerspecs[layer].mua == 0.0)
		&& (In_Ptr->layerspecs[layer].mus == 0.0))
		/* glass layer. */
		HopInGlass(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
	else
		if (In_Ptr->objLayer == 0) /* Conventional boundary check and hop drop spin if the object doesn't exist*/
			HopDropSpinInTissue(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);
		else /** Boundary check involves object boundaries. Hop drop spin depends on Photon_Ptr->inObj **/
			HopDropSpinInTissueObj(In_Ptr, Photon_Ptr,Stokes_Ptr, Out_Ptr);


		Roulette(Photon_Ptr, In_Ptr);

}

void rotSphi(double* S, double phi, double* S2) {
	double	cos2phi, sin2phi;

	cos2phi = cos(2*phi);
	sin2phi = sin(2*phi);

	S2[0] = S[0];
	S2[1] = S[1]*cos2phi+S[2]*sin2phi;
	S2[2] = -S[1]*sin2phi+S[2]*cos2phi;
	S2[3] = S[3];

}
