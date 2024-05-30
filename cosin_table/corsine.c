#include<stdio.h>
#include <math.h>
#include <time.h>
#define ITERS 20
double theta_table[ITERS];


double compute_K(int n) {
    double k = 1.0;
    for (int i = 0; i < n; i++) {
        k *= 1.0 / sqrt(1.0 + pow(2.0, -2.0 * i));
    //printf("k=%.12lf\n",k);
    }
    printf("k=%.12lf\n",k);
    return k;
}


double  CORDIC(double alpha, int n) {
    double K_n = compute_K(n);
    double theta = 0.0;
    double x = 1.0;
    double y = 0.0;
    double P2i = 1.0; // This will be 2**(-i) in the loop below
    for (int i = 0; i < n; i++) {
        double sigma = (theta < alpha) ? 1.0 : -1.0;
        theta += sigma * theta_table[i];
   // printf("theta=%.12lf\n",theta);
	double temp=x;
        x = x - sigma * y * P2i;
        y = sigma * P2i * temp + y;
        P2i /= 2.0;
//	printf("%lf %lf\n",*x,*y);
    }
    return y * K_n;
}

double radians(int x)
{
return x*(M_PI/180);
}
int main() {
	clock_t start,end;
	double cpu_time_used;
	start=clock();
    // Compute the theta_table
    for (int i = 0; i < ITERS; i++) {
        theta_table[i] = atan2(1.0, pow(2.0, i));
	printf("%.8lf, ",theta_table[i]);
    }


    printf("  x       sin(x)     diff. sine    \n");
    //for (int x = -90; x <= 90; x += 15) {
        double sin_x;
	int x=75;
        sin_x=CORDIC(radians(x), ITERS);
        printf("%+05.1f°  %+.8f (%+.8f) \n",
               (double)x, sin_x, sin_x - sin(radians(x)));
    //}
end= clock();
cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC; // in seconds 
 
    printf("main() took %0.12lf seconds to execute \n", cpu_time_used); 
    return 0;
}


