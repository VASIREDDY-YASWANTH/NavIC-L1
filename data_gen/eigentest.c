#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"matfun.h"
#include"L1.h"
int main() {
// r0 register initial parameters for data signal
double **sv_l1_d_r0= createMat(1,64);	sv_l1_d_r0=loadMat("sv_l1_d_r0.dat",1,64); 
// r1 register initial parameters for data signal
double **sv_l1_d_r1= createMat(1,64);	sv_l1_d_r1=loadMat("sv_l1_d_r1.dat",1,64); 
// C register initial parameters for data signals
double **sv_l1_d_c= createMat(1,64);	sv_l1_d_c=loadMat("sv_l1_d_c.dat",1,64); 
// r0 register  initial parameters for pilot signal
double **sv_l1_p_r0= createMat(1,64);	sv_l1_p_r0=loadMat("sv_l1_p_r0.dat",1,64); 
// r1 register initial parameters for pilot signal
double **sv_l1_p_r1= createMat(1,64);	sv_l1_p_r1=loadMat("sv_l1_p_r1.dat",1,64); 
// C register initial patameter for pilot signal
double **sv_l1_p_c= createMat(1,64);	sv_l1_p_c=loadMat("sv_l1_p_c.dat",1,64); 
// r0 register initial parameter for pilot overlay signal
double **sv_l1_ol_r0= createMat(1,64);	sv_l1_ol_r0=loadMat("sv_l1_ol_r0.dat",1,64); 
// r1 register initial parameter for pilot overlay signal
double **sv_l1_ol_r1= createMat(1,64);	sv_l1_ol_r1=loadMat("sv_l1_ol_r1.dat",1,64); 


long int bin = convert(floor(sv_l1_d_r0[0][0]));
printf("\n%ld\n",bin);

// Printing matrices . 	
printf("SV_L1_D_r0 = \n"); 	printMat(sv_l1_d_r0,1,64);
printf("SV_L1_D_r1 = \n");	printMat(sv_l1_d_r1,1,64);
printf("SV_L1_D_C = \n");	printMat(sv_l1_d_c,1,64);
printf("SV_L1_P_r0 = \n");	printMat(sv_l1_p_r0,1,64);
printf("SV_L1_P_r1 = \n");	printMat(sv_l1_p_r1,1,64);
printf("SV_L1_P_c = \n");	printMat(sv_l1_p_c,1,64);
printf("SV_L1_OL_r0 = \n");	printMat(sv_l1_ol_r0,1,64);
printf("SV_L1_OL_r1 = \n");	printMat(sv_l1_ol_r1,1,64);
}

