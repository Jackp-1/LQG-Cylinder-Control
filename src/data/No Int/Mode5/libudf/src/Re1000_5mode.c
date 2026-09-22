#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 5//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738844447013,-0.0112353552989863,0.000277295150341886,-0.000289126600489287,0.000409151545599252,
0.00587402065970237,1.00018623732696,7.71768581542265e-08,-0.000330820546595867,0.000860578637570655,
-0.000145754404832502,-2.9428643517803e-05,0.999723866037778,-0.0174799229323593,-0.000190879670780304,
0.000148701117721479,0.000254826984006231,0.015054969825043,1.00000575404391,1.74881306318324e-05,
2.73226369298749e-05,1.97532872589041e-05,9.70773298781136e-05,5.61464861031742e-05,0.999721201408909
};

static double B[DIM] = {-0.000554648939421054,-0.000411574631716908,-0.000868645450154549,0.000826706133420439,0.00122015986386304};

static double C[DIM] = {0.000337019052925793,-0.000142659000904831,-3.83926410607713e-06,2.51766983532156e-06,2.01567875748752e-05};

static double K[DIM+1] = {-0.0303324122941912,-0.069564502710853,0.00179057970372848,0.000850485787149747,-0.00195510038998623,8.8990565175703e-16};

static double L[DIM] = {61.9863792416346,-59.4073877491547,75.3426544181614,-22.573526567775,-34.745158542865};

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
	fp = fopen("5Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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