#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 6//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738832798355,-0.0112353734367439,0.000277271468965309,-0.000289154072013131,0.000409036053476595,0.000230065434100087,
0.00587400931592537,1.00018621966392,5.41152937461042e-08,-0.000330847299108018,0.00086046616822083,0.000269514796494292,
-0.000145574842782777,-2.91490530967102e-05,0.999724231082131,-0.0174794994635325,-0.000189099379522089,-0.00230355059762498,
0.000148668843571372,0.000254776730941915,0.0150549042126473,1.00000567793041,1.71681443792678e-05,0.000618138306489641,
2.83107881547907e-05,2.12919064576326e-05,9.9086212186056e-05,5.84768852075733e-05,0.999730998562543,-0.0127338926767205,
-7.56444865894026e-05,-0.000118678912466618,0.000267373513182962,-0.000256974560446173,0.0102216270507081,0.999775644563646
};

static double B[DIM] = {-0.000554893786217863,-0.000411813070144719,-0.000864871178774902,0.000826027752875758,0.00124093012624746,-0.00152888126966396};

static double C[DIM] = {0.000337020976120357,-0.00014265853006757,-3.83465358122336e-06,2.51365026633519e-06,2.01519462415884e-05,-3.88848662207161e-05};

static double K[DIM+1] = {-0.114871184614992,-0.0394625081445206,0.0027399745569569,0.000692961002811338,0.00806076242810546,0.0217876659750245,-1.56170348322076e-31};

static double L[DIM] = {57.9389564785523,-57.144626768589,49.4442453756649,-17.8514998243037,-28.5552902761438,-50.8231016728122};

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
	fp = fopen("6Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], xhat[6], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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