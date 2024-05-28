#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 131072

typedef struct {
    double real;
    double imag;
} complex;

int Log2(int N) {
    int k = N, i = 0;
    while (k) {
        k >>= 1;
        i++;
    }
    return i - 1;
}

int check(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n) {
    int j, p = 0;
    for (j = 1; j <= Log2(N); j++) {
        if (n & (1 << (Log2(N) - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

void ordina(complex* f1, int N) {
    complex f2[MAX];
    for (int i = 0; i < N; i++)
        f2[i] = f1[reverse(N, i)];
    for (int j = 0; j < N; j++)
        f1[j] = f2[j];
}

void transform(complex* f, int N) {
    ordina(f, N);
     //   for(int i=0;i<N;i++)
	//	printf("%lf + %lfi\n", f[i].real, f[i].imag);



    complex* W;
    W = (complex*)malloc(N / 2 * sizeof(complex));
    W[1].real = cos(-2. * M_PI / N);
    W[1].imag = sin(-2. * M_PI / N);
    W[0].real = 1;
    W[0].imag = 0;
	//	printf("%lf + %lfi\n", W[1].real, W[1].imag);
    
    for (int i = 2; i < N / 2; i++) {
        double angle = -2. * M_PI * i / N;
        W[i].real = cos(angle);
        W[i].imag = sin(angle);
//	printf("%lf , %lf\n",W[i].real,W[i].imag);
    }
    printf("\n");

    int n = 1;
    int a = N / 2;
       // for(int i=0;i<N;i++)
//		printf("%lf + %lfi\n", f[i].real, f[i].imag);
//	printf("\n");

    for (int j = 0; j < log2(N); j++) {
        for (int i = 0; i < N; i++) {
            if (!(i & n)) {
                complex temp = f[i];
                complex Temp;
                Temp.real = W[(i * a) % (n * a)].real * f[i + n].real - W[(i * a) % (n * a)].imag * f[i + n].imag;
                Temp.imag = W[(i * a) % (n * a)].real * f[i + n].imag + W[(i * a) % (n * a)].imag * f[i + n].real;
		//printf("%lf + %lfi\n", temp.real, temp.imag);
		//printf("%lf + %lfi\n\n", Temp.real, Temp.imag);
                f[i].real = temp.real + Temp.real;
                f[i].imag = temp.imag + Temp.imag;
                f[i + n].real = temp.real - Temp.real;
                f[i + n].imag = temp.imag - Temp.imag;
            }
        }
        n *= 2;
        a = a / 2;
    }
    free(W);
}

void FFT(complex* f, int N, double d) {
    transform(f, N);
    for (int i = 0; i < N; i++) {
        f[i].real *= d;
        f[i].imag *= d;
    }
}

int main() {
    int a;
        scanf("%d", &a);
    //int n=131072;
    int n=pow(2,a);
    double d=1;
    FILE *file;
    complex vec[MAX];

file=fopen("dummy.dat","r");
//for(int i=0;i<n;i++)
for(int i=0;i<131072;i++)
{
	vec[i].imag=0;
	fscanf(file ,"%lf",&vec[i].real);
}
fclose(file);		

//for (int i = 102300; i < 131072; i++) {
//	vec[i].imag=0;
//	vec[i].real=0;
//}

    FFT(vec, n, d);
    printf("...printing the FFT of the array specified\n");
    for (int j = 0; j < n; j++)
    {
        printf("%lf + %lfi\n", vec[j].real, vec[j].imag);
    }


FILE *fp=fopen("fftout1.dat","w");
if(fp==NULL)
	printf("error file not found");
for (int i = 0; i < n; i++) {
    fprintf(fp, "%.3lf + %.3lfi\n", vec[i].real, vec[i].imag);
  }
fclose(fp);

    

return 0;
}


