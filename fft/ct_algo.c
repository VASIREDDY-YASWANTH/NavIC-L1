#include<stdio.h>
#include<stdbool.h>
#include<complex.h>
#include<math.h>
bool check(int n)
{
if(n==0)
	return false;
else 
	return (ceil(log2(n))==floor(log2(n)));
}







int main()
{

long int n;  // array dimension
do{ printf("enter array dimension 'n' (must be a power of 2)");
    scanf("%ld",&n);
}while(!check(n));

int d=1;   // sampling step size

double complex arr[200];
printf("enter array values\n");
for(int i=0;i<n;i++)
{
	double re=0,im=0;// real and imaginary parts
	printf("enter the array element %d =", i);
	scanf("%lf", &re  ); // reading only real part taking img part always zero 
	arr[i]=CMPLX(re,im);
}

for(int i=0;i<n;i++)
printf("%.2lf +%.0lf i \n",creal(arr[i]) , cimag(arr[i]));






}
