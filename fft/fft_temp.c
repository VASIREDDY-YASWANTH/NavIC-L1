
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 200

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
    complex* W;
    W = (complex*)malloc(N / 2 * sizeof(complex));
    W[1].real = cos(-2. * M_PI / N);
    W[1].imag = sin(-2. * M_PI / N);
    W[0].real = 1;
    W[0].imag = 0;
    for (int i = 2; i < N / 2; i++) {
        double angle = -2. * M_PI * i / N;
        W[i].real = cos(angle);
        W[i].imag = sin(angle);
    }
    int n = 1;
    int a = N / 2;
    for (int j = 0; j < log2(N); j++) {
        for (int i = 0; i < N; i++) {
            if (!(i & n)) {
                complex temp = f[i];
                complex Temp;
                Temp.real = W[(i * a) % (n * a)].real * f[i + n].real - W[(i * a) % (n * a)].imag * f[i + n].imag;
                Temp.imag = W[(i * a) % (n * a)].real * f[i + n].imag + W[(i * a) % (n * a)].imag * f[i + n].real;
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
    int n;
    FILE *file;
    
file=fopen("subcarrier.dat","r");
for(i=0;i<siglen;i++)
	fscanf(file ,"%f",&subCarrier[i]);
fclose(file);		
    do {
        printf("specify array dimension (MUST be power of 2)\n");
        scanf("%d", &n);
    } while (!check(n));
    double d;
    printf("specify sampling step\n");
    scanf("%lf", &d);
    complex vec[MAX];

   // printf("specify the array\n");

    //for (int i = 0; i < n; i++) {
     //   printf("specify element number: %d\n", i);
	//vec[i].imag=0;
      //  scanf("%lf", &vec[i].real);
    //}

    FFT(vec, n, d);
    printf("...printing the FFT of the array specified\n");
    for (int j = 0; j < n; j++)
    {
        printf("%lf + %lfi\n", vec[j].real, vec[j].imag);
    }
    return 0;
}


