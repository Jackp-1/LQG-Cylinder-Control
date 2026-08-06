#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 4//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738828571859,-0.0112353791688803,0.000277261029267997,-0.000289162187084437,
0.0058739879469176,1.00018618813999,6.86602042235251e-09,-0.000330893877327675,
-0.000145750236900131,-2.94223766112503e-05,0.999723874996074,-0.0174799135892986,
0.000148697526116201,0.00025482158366559,0.0150549621054686,1.00000574599277
};

static double B[DIM] = {-0.000555190443883629,-0.000412690470934233,-0.000868503281200388,0.000826583623092767};

static double C[DIM] = {0.000337018560960106,-0.000142659121348067,-3.84044350834572e-06,2.51869806723247e-06};

static double K[DIM+1] = {-0.0389250091351826,-0.0551709564755442,0.00179079461972592,0.000765058912112125,-0};

static double L[DIM] = {60.086378759469,-56.1987330430438,73.5626369411958,-25.5281755498168};

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
	fp = fopen("4Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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