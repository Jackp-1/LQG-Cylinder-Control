#include "udf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID_cyl 8
#define Ts 0.4
#define RAD 0.5
#define ref 0
#define DIM 10//TRIPLE CHECK THIS IS CORRECT

static double cyl;

static double A[DIM][DIM] = {0.999738751147911,-0.0112352944340849,0.000277050135968115,-0.000288509895127518,0.00041084823141434,0.000221570146210398,-0.00120361225998704,0.000451782567723562,0.00162840791843139,-0.000964795983335201,
0.00587401199018992,1.00018621770063,5.49073883894524e-08,-0.00033087230316299,0.000860425420275398,0.000269526986654137,-9.22585119655425e-05,-6.53004389184714e-05,-0.00012208153477721,0.000303281988786478,
-0.000145396681727377,-2.9210771422733e-05,0.999724766643997,-0.0174803541213381,-0.000191980930678638,-0.00228767518475992,1.72111752519682e-06,-0.00443117184122434,0.000629337525615888,7.42453990909848e-05,
0.00014878052948248,0.000254819794497574,0.0150551887758729,1.00000544093939,1.62785817879125e-05,0.000623737949494867,0.000699489291035245,-0.0015967900499476,0.000723610362786736,0.000877879040729343,
2.94286539480135e-05,2.18435220410715e-05,0.000102070929448446,5.67066294555666e-05,0.999722880957654,-0.0126725333619375,0.024453087122601,0.00400343503437767,-0.00095365008785189,-0.00360262666904872,
-7.68139576875228e-05,-0.00011882082420603,0.000264000940773688,-0.000253589044136263,0.0102347248479629,0.999694332320544,-0.00602716848821133,0.0287582621475463,-0.00940560482411887,0.000451104835072693,
7.25639740374289e-05,2.23720330957459e-05,0.000135532549942031,-0.000219880853215887,-0.0182670072810978,0.00563837871387232,0.999422250168191,0.0144202363674049,-0.00295935281277219,-0.000363163887211842,
-4.09571705385972e-05,2.67287418213285e-05,0.00052045161810468,0.00044546542078554,-0.00260883315057279,-0.0253538358488434,-0.0130922026218363,0.999240380482307,-0.00508177087590298,-0.00327445941504363,
-8.8491788795015e-06,6.19584172959139e-05,-6.40295731816394e-05,8.68511706765429e-05,0.000825759769447438,0.00645376730516347,0.00291612313589627,0.00344370704892918,0.997484918919285,0.003767449552558,
2.59166046791531e-05,7.09949113033952e-06,2.59478071550127e-05,-0.000114245282883471,0.000893338463710079,-2.81726049143386e-05,-0.000387742484150536,0.0023525675280949,0.00221307640116484,0.99646358025902
};

static double B[DIM] = {-0.000704815010901108,-0.000393840100664469,-0.000825841565917376,0.000855446763851523,0.00116980585302589,-0.00135723721575536,0.00730746896597744,-0.0151220446431838,-0.0484831859598124,0.0547591409352706};

static double C[DIM] = {0.000337023994908255,-0.000142657791006711,-3.82741656168501e-06,2.50734085467095e-06,2.01443469281262e-05,-3.88659940626897e-05,3.58286334101238e-05,4.69256337128748e-05,9.74213031987081e-05,-8.72339109186023e-05};

static double K[DIM+1] = {-0.0319431882846607,0.0753998546573229,-0.000860992788833238,-0.000581165810513397,0.00188432642210782,0.00123662392870708,0.00196153801896393,-0.00253830875743502,-0.0137692259374524,0.00817968123794331,2.21203309424936e-14};

static double L[DIM] = {59.5582228865686,-53.4670825783328,73.1883601091543,38.0312248065769,-49.4903404574768,-103.941502078192,12.3459802753895,15.1287413556576,18.3564709526352,8.16495878253907};

static double Unew = 0.0; //1e-5 end data

static double U = 0; 

static double xhat[DIM+1] = {0.0};

static double integral = 0.; // 1e-5 end integral data


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

Unew = -Kxhat;

if I_AM_NODE_ZERO_P {
	fp = fopen("10Mode.txt","a");
	fprintf(fp,"%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n", xhat[0] ,xhat[1] ,xhat[2], xhat[3], xhat[4], xhat[5], xhat[6], xhat[7], xhat[8], xhat[9], xhat[10], y , yhat, noise_dub, y_noise, Unew, force_xcyl);
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

	F_PROFILE(f, thread, position) = 0.01460734693878+((0.01460734693878*0.05)*sin((0.0203/5)*t));
	}
    end_f_loop(f, thread)
}