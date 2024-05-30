#include<stdio.h>
#include <math.h>

#define ITERS 32

double theta_table[ITERS];

double compute_K(int n) {
    double k = 1.0;
    for (int i = 0; i < n; i++) {
        k *= 1.0 / sqrt(1.0 + pow(2.0, -2.0 * i));
    }
   // printf("k=%.12lf\n",k);
    return k;
}

void CORDIC(double alpha, int n, double *x, double *y) {
    double K_n = compute_K(n);
    double theta = 0.0;
    *x = 1.0;
    *y = 0.0;
    double P2i = 1.0; // This will be 2**(-i) in the loop below
    for (int i = 0; i < n; i++) {
        double sigma = (theta < alpha) ? 1.0 : -1.0;
        theta += sigma * theta_table[i];
	double temp=*x;
        *x = *x - sigma * *y * P2i;
        *y = sigma * P2i * temp + *y;
        P2i /= 2.0;
//	printf("%lf %lf\n",*x,*y);
    }
    *x *= K_n;
    *y *= K_n;
}

double radians(int x)
{
return x*(M_PI/180);
}
int main() {
    // Compute the theta_table
    for (int i = 0; i < ITERS; i++) {
        theta_table[i] = atan2(1.0, pow(2.0, i));
	printf("%lf\n",theta_table[i]);
    }


    printf("  x       sin(x)     diff. sine     cos(x)    diff. cosine \n");
    for (int x = -90; x <= 90; x += 15) {
        double cos_x, sin_x;
        CORDIC(radians(x), ITERS, &cos_x, &sin_x);
        printf("%+05.1f°  %+.8f (%+.8f) %+.8f (%+.8f)\n",
               (double)x, sin_x, sin_x - sin(radians(x)), cos_x, cos_x - cos(radians(x)));
    }

    return 0;
}


