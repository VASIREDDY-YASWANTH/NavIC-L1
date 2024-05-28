#include<stdio.h>
#include<stdlib.h>
float readfile(char *filename, int arrsize);
int main()
{
int siglen=102300;
float subCarrier[siglen];
float pilotSig[siglen];
float pilotXsubc[siglen];
int i=0;
FILE *file;
// reading carrier signal 
file=fopen("subcarrier.dat","r");
for(i=0;i<siglen;i++)
	fscanf(file ,"%f",&subCarrier[i]);
fclose(file);		
	printf("carriersig=");
	for(i=0;i<20;i++)
	printf("%0.1f  " , subCarrier[i]);
	printf("\n");

//reaading pilot signal
file=fopen("pilotsig.dat","r");
for(i=0;i<siglen;i++)
	fscanf(file ,"%f",&pilotSig[i]);
fclose(file);		
	printf("pilotsig  =");
	for(i=0;i<20;i++)
	printf("%0.1f  " , pilotSig[i]);
	printf("\n");
// pilot signal X sub carrier
for (int i=0;i<siglen;i++)
pilotXsubc[i]=subCarrier[i]*pilotSig[i];
	printf("pilotXsubc=");
	for(i=0;i<20;i++)
	printf("%0.1f  " , pilotSig[i]);
	printf("\n");


}









/*
float readfile(char *filename ,int arrsize)
{

	float arr[arrsize];
	FILE *file;
	file =  fopen(filename ,"r");
	for (int i=0;i<arrsize;i++)
	{
		fscanf(file,"%f",&arr[i]);

	}
	fclose(file);
// return arr;
}

*/
