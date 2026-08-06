#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 9//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738771363071,-0.0112352843802957,0.000277078794763826,-0.000288567574398714,0.000410750112631646,0.000221717900910865,-0.00120415940766643,0.000452980822922064,0.00163302891213259,
0.00587400593880196,1.00018621469104,4.63284063126346e-08,-0.000330855036930345,0.000860454792035669,0.000269482756431158,-9.20947238504974e-05,-6.56591354234589e-05,-0.000123464824657553,
-0.000145397859017073,-2.92113569349099e-05,0.999724764974967,-0.0174803507622149,-0.000191975216440861,-0.00228768378969261,1.75298229093711e-06,-0.00443124162516511,0.000629068408360203,
0.000148762532521619,0.000254810843905301,0.0150551632617915,1.00000549228955,1.63659340491121e-05,0.000623606407838644,0.00069997640049268,-0.00159785682124453,0.000719496428250959,
2.95045310414084e-05,2.18812586861789e-05,0.000102178499518141,5.64901317601271e-05,0.999722512671268,-0.0126719787683814,0.0244510334175321,0.00400793265516226,-0.000936305304084383,
-7.68247993335842e-05,-0.000118826216180489,0.000263985570698588,-0.00025355811001264,0.0102347774703086,0.999694253077827,-0.00602687504598171,0.0287576195081084,-0.00940808312156187,
7.25747543294563e-05,2.23773945564254e-05,0.000135547833036147,-0.000219911612279841,-0.0182670596056479,0.00563845750814607,0.999421958386583,0.0144208753700757,-0.00295688854028617,
-4.08971530949359e-05,2.67585908425684e-05,0.000520536704134077,0.000445294174923847,-0.00260912445861358,-0.0253533971751143,-0.0130938270663529,0.999243938021085,-0.0050680514576144,
-8.97082419254306e-06,6.18979183257856e-05,-6.42020283221392e-05,8.71982573762415e-05,0.000826350202089026,0.00645287818693669,0.00291941561304313,0.00343649651326219,0.997457111954674
};

static double B[DIM] = {-0.000651572209417741,-0.000409778280470868,-0.000828942318168598,0.000808046266476655,0.0013696513708115,-0.00138579200446634,0.00733586216000806,-0.0149639703605459,-0.0488035760733606};

static double C[DIM] = {0.000337023322521602,-0.00014265795562068,-3.82902849187389e-06,2.5087461750632e-06,2.01460395534839e-05,-3.88701975336201e-05,3.58334200082917e-05,4.69221630715552e-05,9.74238404383842e-05};

static double K[DIM+1] = {-0.0157627865968722,0.0845498558072999,-0.00133009066911557,-0.000707021300832762,0.00178941692135074,-0.00010139788711797,0.00382932484533442,-0.00217475076243921,-0.0157439422258635,-0};

static double L[DIM] = {60.2813368706063,-54.4408904960757,79.8416104908371,25.3544126180031,-45.7739221445924,-92.5985148017697,9.74224269692398,28.4380596970432,-1.8389685389395};

static double Unew = 0.0; //1e-5 end data

static double U = 0; 

static double xhat[DIM+1] = {0.0};

static double time = 0.;

static double rate;

static double integral = 0.; // 1e-5 end integral data

static double iter = 0.;

static int count = 0;

static double rate_u = 3e-4;

static double rate_d = -3e-4;

static double limit = 2*0.01460734693878;

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
double y = force_ycyl; // Measurement from Simulation
double yhat = 0;
double y_noise = y;//+noise_dub;

U = Unew;
  
for (i = 0; i < DIM; i++) {
    yhat += C[i] * xhat[i];
}
integral += (yhat - ref);

for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
        Axhat[i] += A[i][j] * xhat[j];
    }
}

for (i = 0; i < DIM; i++) {
    Bu[i] = B[i] * U;
}

for (i = 0; i < DIM; i++) {
	xhat[i] = Axhat[i] + Bu[i] + L[i] * (y_noise - yhat);
}

xhat[DIM] = -integral;
	
for (i = 0; i < DIM+1; i++) {
    Kxhat += K[i] * xhat[i];
}

//if (iter<300.0){
//   Unew = 0;
//}
//else{
	Unew = -Kxhat;
//	if (Unew>limit){
//		Unew = limit;
//	} else if (Unew<-limit){
//		Unew = -limit;
//	}else{
//			Unew = -Kxhat;
//	}
//}

iter++;

if I_AM_NODE_ZERO_P {
	fp = fopen("9Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], xhat[6], xhat[7], xhat[8], xhat[9], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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