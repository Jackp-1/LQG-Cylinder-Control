#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 7//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738766483056,-0.0112353975567217,0.000277095568032643,-0.000289011298113822,0.000409600517455712,0.000226242316209163,-0.00120932663916957,
0.00587400513410513,1.00018621814293,4.30230407619248e-08,-0.000330838295834712,0.000860501763127457,0.000269273711978635,-9.17782261601382e-05,
-0.000145564938088829,-2.91454505947086e-05,0.999724257354273,-0.0174795207878995,-0.00018918368650732,-0.00230297958591219,3.30775490480602e-05,
0.000148703117886995,0.000254789197080846,0.0150549951250635,1.00000560413932,1.68764075289972e-05,0.000620114241916878,0.000708733367338724,
2.96819828138374e-05,2.17906327709414e-05,0.000102723297779244,5.55247639008019e-05,0.999719327198604,-0.0126548424592444,0.0244166741254057,
-7.58124758619136e-05,-0.000118740012960269,0.000266927922647159,-0.000256612886990507,0.0102230569453802,0.999765959878515,-0.00619101294378862,
7.33506829664115e-05,2.19898387315396e-05,0.000137927634134878,-0.00022409735731925,-0.0182808952245144,0.00571301155953393,0.99927217896125
};

static double B[DIM] = {-0.000563668013319825,-0.000412366369925579,-0.00086356068127238,0.000830562613341102,0.00142235392232196,-0.00155110805658818,0.00755957670366995};

static double C[DIM] = {0.000337022810954806,-0.000142658080862667,-3.83025488439199e-06,2.5098153742297e-06,2.01473273406901e-05,-3.88733956282612e-05,3.58370617591691e-05};

static double K[DIM+1] = {-0.0242228013415966,-0.0556441622110824,0.0015157038116844,0.000783278248563715,-0.00208732319467441,-0.00128357484940455,-0.00219017926156968,9.60673050688754e-16};

static double L[DIM] = {61.2812603284712,-57.0860908413589,54.7939834841317,-20.218404568921,-53.8286830261913,-43.1939020749326,37.5225295772579};

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
	fp = fopen("7Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], xhat[6], xhat[7], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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