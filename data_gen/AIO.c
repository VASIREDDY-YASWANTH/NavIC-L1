#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include"NavicL1.h"
int main()
{
int r = 1 , c = 64 ,msl=19; // r=rows,c=columns,msl=maximum string length in array

//************loading arrays from .dat files**************
// r0 register initial parameters for data signal
    char ***sv_l1_d_r0 = loadStrArr("sv_l1_d_r0.dat", r, c, 19);
// r1 register initial parameters for data signal
    char ***sv_l1_d_r1 = loadStrArr("sv_l1_d_r1.dat", r, c, msl);
// C register initial parameters for data signals
    char ***sv_l1_d_c = loadStrArr("sv_l1_d_c.dat", r, c, msl);
// r0 register  initial parameters for pilot signal
    char ***sv_l1_p_r0 = loadStrArr("sv_l1_p_r0.dat", r, c, msl);
// r1 register initial parameters for pilot signal
    char ***sv_l1_p_r1 = loadStrArr("sv_l1_p_r1.dat", r, c, msl);
// C register initial patameter for pilot signal
    char ***sv_l1_p_c = loadStrArr("sv_l1_p_c.dat", r, c, msl);
// r0 register initial parameter for pilot overlay signal
    char ***sv_l1_ol_r0 = loadStrArr("sv_l1_ol_r0.dat", r, c, msl);
// r1 register initial parameter for pilot overlay signal
    char ***sv_l1_ol_r1 = loadStrArr("sv_l1_ol_r1.dat", r, c, msl);

//*************Octal To Binary Conversion of arrays********
char ***r0_data=OctToBin(sv_l1_d_r0,r,c,msl);
char ***r1_data=OctToBin(sv_l1_d_r1,r,c,msl);
char ***r0_pilot=OctToBin(sv_l1_p_r0,r,c,msl);
char ***r1_pilot=OctToBin(sv_l1_p_r1,r,c,msl);


int satId[] = {1, 46, 62, 50};
int satIdlength = sizeof(satId) / sizeof(satId[0]);
genNavicCaTable_Data(4e6, 10230, 1.023e6, satId, satIdlength,sv_l1_d_r0,sv_l1_d_r1,sv_l1_d_c);
genNavicCaTable_Pilot(4e6, 10230, 1.023e6, satId, satIdlength,sv_l1_p_r0,sv_l1_p_r1,sv_l1_p_c);
genNavicCaTable_Pilot_Overlay(4e6, 10230, 1.023e6, satId, satIdlength,sv_l1_ol_r0,sv_l1_ol_r1);



//	printf("%s\n",r0_data[0][0]);


//****************printing arrays*********************
//printf("SV_L1_D_r0 = \n");	printStrArr(sv_l1_d_r0, r, c);
//printf("SV_L1_D_r1 = \n");	printStrArr(sv_l1_d_r1, r, c);
//printf("SV_L1_D_c = \n");	printStrArr(sv_l1_d_c, r, c);
//printf("SV_L1_P_r0 = \n");	printStrArr(sv_l1_p_r0, r, c);
//printf("SV_L1_P_r1 = \n");	printStrArr(sv_l1_p_r1, r, c);
//printf("SV_L1_P_c = \n");	printStrArr(sv_l1_p_c, r, c);
//printf("SV_L1_OL_r0 = \n");	printStrArr(sv_l1_ol_r0, r, c);
//printf("SV_L1_OL_r1 = \n");	printStrArr(sv_l1_ol_r1, r, c);
//printf("r0_data = \n");	printStrArr(r0_data, r, c);
//printf("r1_data = \n");	printStrArr(r1_data, r, c);
//printf("r0_pilot = \n");	printStrArr(r0_pilot, r, c);
//printf("r1_pilot = \n");	printStrArr(r1_pilot, r, c);
return 0;
}


//printing prn data of sattelite ids
/*
  	for(int i=0;i<64;i++)
  	{
	printf("\nSATELLITE: %d\n",i+1);
  	//int* data_prn=genNavicCaCode_Data(i,sv_l1_d_r0,sv_l1_d_r1,sv_l1_d_c);
  	int* data_prn=genNavicCaCode_Data(i,r0_data,r1_data,sv_l1_d_c);
  	printf("Data prn (first 24)          :");
  	for(int j=0;j<24;j++){
  		printf("%d",data_prn[j]);
  	}
  	
  	printf(" Data prn (Last 24)           :");
  	for(int j=10206;j<10230;j++){
  		printf("%d",data_prn[j]);
  	}
  	//int* pilot_prn=genNavicCaCode_Pilot(i,sv_l1_p_r0,sv_l1_p_r1,sv_l1_p_c);
  	int* pilot_prn=genNavicCaCode_Pilot(i,r0_pilot,r1_pilot,sv_l1_p_c);
  	printf("\nPilot prn (first 24)         :");
  	for(int j=0;j<24;j++){
  		printf("%d",pilot_prn[j]);
  	}
  	printf(" Pilot prn (last 24)          :");
  	for(int j=10206;j<10230;j++){
  		printf("%d",pilot_prn[j]);
  	}
  	
  	//int* pilot_overlay_prn=genNavicCaCode_Pilot_Overlay(i,sv_l1_ol_r0,sv_l1_ol_r1);
  	int* pilot_overlay_prn=genNavicCaCode_Pilot_Overlay(i,sv_l1_ol_r0,sv_l1_ol_r1);
  	printf("\npilot overlay prn (first 24) :");
  	for(int j=0;j<24;j++){
  		printf("%d",pilot_overlay_prn[j]);
  	}
  	
  	printf(" pilot overlay prn (last 24)  :");
  	for(int j=1776;j<1800;j++){
  		printf("%d",pilot_overlay_prn[j]);
  	}
  }
*/
