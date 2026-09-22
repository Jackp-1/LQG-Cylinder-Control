#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 8//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738751670448,-0.0112353873901585,0.00027704937728935,-0.000288920137500518,0.00040989424901676,0.000224699515771562,-0.00120632105579246,0.000457765043053614,
0.00587400642882597,1.0001862172543,4.70604198839208e-08,-0.000330846263880449,0.000860476089029562,0.000269408563028815,-9.20409341973569e-05,-6.57781841891668e-05,
-0.000145406100638915,-2.92544679178497e-05,0.999724752663367,-0.0174804983145119,-0.000192333406624645,-0.00228643594464263,8.4830410819795e-07,-0.00442923936609142,
0.000148751982692817,0.000254755658955429,0.0150551475021175,1.00000530341271,1.59074265874306e-05,0.000625203733204158,0.000698818351735746,-0.00159529379532495,
2.95097982447494e-05,2.19088108253966e-05,0.000102186367836011,5.65844321251511e-05,0.999722741589867,-0.0126727762635483,0.0244516115954705,0.00400665301556725,
-7.66928108099016e-05,-0.000118135799284392,0.000264182739405304,-0.000251195078785582,0.0102405138404442,0.999674268998418,-0.00601238674020723,0.0287255535842776,
7.25992588787774e-05,2.25055750557273e-05,0.000135584438720937,-0.000219472899779859,-0.0182659946102919,0.005634747330584,0.999424648237133,0.0144149221175553,
-4.08127043810798e-05,2.72003324219214e-05,0.000520662856342106,0.000446806086280955,-0.00260545422224206,-0.0253661833610405,-0.0130845571786785,0.999223421646244
};

static double B[DIM] = {-0.000570077073420304,-0.000411806175318315,-0.00079483553147701,0.00085170524121854,0.00134785379462258,-0.00193200784562799,0.007234453547857,-0.0153134494212641};

static double C[DIM] = {0.00033702116216182,-0.000142658484520824,-3.83420757913551e-06,2.51326143074328e-06,2.01514779120982e-05,-3.88837031698442e-05,3.58487992162545e-05,4.69110119966691e-05};

static double K[DIM+1] = {-0.0265036690030147,-0.0815719271424728,0.00194802331074451,0.000996284716048743,-0.0012812622399229,0.0017492345853954,-0.000772817762809147,-0.00264604293174599,-0};

static double L[DIM] = {62.6750763343467,-55.0818455534548,77.732355515608,-8.37412021924097,-44.7697029329625,-90.5513244551211,11.4370919200376,28.3011380411946};

static double Unew = 0;

static double U = 0; 

static double xhat[DIM+1] = {0.0};

static double time = 0;

static double rate;

static double integral = 0.0;

static double iter = 0;

static int count = 0;

static double rate_u = 3e-4;

static double rate_d = -3e-4;

static double limit = 0.0;

static double iter2 = 0.0;

DEFINE_EXECUTE_AT_END(Control_loop)
{
FILE *fp;

Domain *d = Get_Domain(1); // for single phase flow

Thread *t_objectcyl = Lookup_Thread(d, ID_cyl);

real forcecyl[ND_ND], moment[ND_ND], cg[ND_ND]; //initialise

Compute_Force_And_Moment (d,t_objectcyl,cg,forcecyl,moment,TRUE);

real force_xcyl = forcecyl[0]; //force components of surface "Boundary_ID" 
real force_ycyl = forcecyl[1];

int upper_bound = 10000;
int noise = (rand() %upper_bound)+1;

double noise_dub = ((noise-(upper_bound/2))*5e-4)*2e-5; // multiply by 1/(2*upper_bound) now between -1 and 1 then mult by noise variance scaling

int i, j;

double Axhat[DIM] = { 0 };
double Bu[DIM] = { 0 };
double Kxhat = 0;
double cl = 0;
double y = force_ycyl; // Measurement from Simulation
double yhat = 0;
double y_noise = y;//+noise_dub;

 
  U = Unew;
  
 for (i = 0; i < DIM; i++) {
    yhat += C[i] * xhat[i];
}
integral = integral + (yhat - ref);

for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
        Axhat[i] += A[i][j] * xhat[j];
    }
}

for (i = 0; i < DIM; i++) {
    Bu[i] = B[i] * Unew;
}

for (i = 0; i < DIM; i++) {
	xhat[i] = Axhat[i] + Bu[i] + L[i] * (y_noise - yhat);
}

xhat[DIM] = integral;
	
for (i = 0; i < DIM+1; i++) {
    Kxhat += K[i] * xhat[i];
}

//if (iter<300.0){
//   Unew = 0;

//}

//else{
	Unew = -Kxhat;
//}

iter++;

if I_AM_NODE_ZERO_P {
	fp = fopen("8Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], xhat[6], xhat[7], xhat[8], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
	fclose(fp);
}
}

DEFINE_PROFILE(update_top, thread, position)
{
    face_t f;
    real vel, x[ND_ND];   
// update control velocity
    begin_f_loop(f, thread)
    {
		F_CENTROID(x, f, thread);  	 
		if(Unew>0){
			F_PROFILE(f, thread, position) = fabs(Unew);
		}
		else{
			F_PROFILE(f, thread, position) = 0;
		}
    }
    end_f_loop(f, thread)
}

DEFINE_PROFILE(update_bot, thread, position)
{
    face_t f;
    real vel, x[ND_ND];   
// update control velocity
    begin_f_loop(f, thread)
    {
		F_CENTROID(x, f, thread);  	 
		if(Unew<0){
			F_PROFILE(f, thread, position) = fabs(Unew);
		}
		else{
			F_PROFILE(f, thread, position) = 0;
		}
    }
    end_f_loop(f, thread)
}

DEFINE_PROFILE(inlet_vel, thread, position)
{
    face_t f;
    real vel, x[ND_ND];
    real t = CURRENT_TIME;

// update control velocity
    begin_f_loop(f, thread)
    {
	F_CENTROID(x, f, thread); 
	//F_PROFILE(f, thread, position) = 0.01460734693878;
	//if(iter>900.0){
	F_PROFILE(f, thread, position) = 0.01460734693878+((0.01460734693878*0.05)*sin((0.0203/5)*t));
    //}
	}
    end_f_loop(f, thread)
}