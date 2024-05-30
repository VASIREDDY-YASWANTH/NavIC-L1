#include<stdio.h>
#include<math.h>

int itr=16;

float compute_k(int n)
{
float k=1;
for (int i=0;i<n;i++)
	k*=1/sqrt(1+pow(2,(-2*i)));

return k;
}

float cordic(float alpha ,int n,float theta_table[])
{
float k_n=compute_k(n);
float theta=0,x=1,y=0;
int p2i=1,sigma;
for(int i=0;i<itr;i++)
{	if(theta<alpha)
		sigma =1;
	else sigma =-1;

theta+=sigma*theta_table[i];
	x=x-sigma*y*p2i;
	y=sigma*p2i*x+y;
	p2i/=2;
}
return x*k_n;
}


int main()
{
int i;
float cos_x,sin_x ,ang;

float theta_table[16] ;
for (int j=0;j<itr;j++) 
theta_table[i]=atan2(1,pow(2,i));


for(i=-90;i<91;i+=15)
{
ang=i*(M_PI/180);
cos_x=cordic(ang,itr,theta_table);
printf("%d   %f  %f\n",i,cos_x,cos(ang));
}


}


