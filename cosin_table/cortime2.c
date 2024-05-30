#include<stdio.h>
#include <math.h>
#include <time.h>
#define ITERS 20

int main() {
clock_t start,end;
double cpu_time_used;
double theta_table[ITERS]={  0.78539816, 0.46364761, 0.24497866, 0.12435499, 0.06241881,
			     0.03123983, 0.01562373, 0.00781234, 0.00390623, 0.00195312, 
			     0.00097656, 0.00048828, 0.00024414, 0.00012207, 0.00006104, 
			     0.00003052, 0.00001526, 0.00000763, 0.00000381, 0.00000191};
double K_n=0.607252935009017017;
double sin_x,alpha,theta,x,y,P2i,sigma,temp;
 int n=ITERS,i,a;

start=clock();
printf("  a       sin(a)     diff. sine    \n");
for ( a = -90; a <= 90; a += 15) {
	//a=79;
	alpha= a*(M_PI/180);
	theta = 0.0; x = 1.0; y = 0.0; P2i = 1.0;//initializing to default values 
    for ( i = 0; i < n; i++) {
        sigma = (theta < alpha) ? 1.0 : -1.0;
        theta += sigma * theta_table[i];
        //printf("theta=%.12lf\n",theta);
	temp=x;
        x = x - sigma * y * P2i;
        y = sigma * P2i * temp + y;
        P2i /= 2.0;
    }
        sin_x=y*K_n;
        //printf("%+05.1f°  %+.8f  \n",(double)a, sin_x);
        printf("%+05.1f°  %+.8f (%+.8f) \n",(double)a, sin_x, sin_x - sin(alpha));
}

end= clock();
cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC; // in seconds 
printf("main() took %0.5lf micro seconds to execute \n", cpu_time_used*1000000); 
return 0;
}
