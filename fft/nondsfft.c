#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MAX 131072

int Log2(int N) {
    int k = N, i = 0;
    while (k) {
        k >>= 1;
        i++;
    }
    return i - 1;
}


int reverse(int N, int n) {
    int j, p = 0;
    for (j = 1; j <= Log2(N); j++) {
        if (n & (1 << (Log2(N) - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

void transform(double fr[],double fi[], int N) {
    //ordina function task
    double fr2[MAX],fi2[MAX];
    for (int i = 0; i < N; i++)
    {
        fr2[i] = fr[reverse(N, i)];
        fi2[i] = fi[reverse(N, i)];
    }
    for (int j = 0; j < N; j++)
    {
        fr[j] = fr2[j];
        fi[j] = fi2[j];
    }

        for(int i=0;i<N;i++)
		printf("%lf + %lfi\n", fr[i], fi[i]);

printf("Log2%d\n",Log2(16));

//for(int i=0;i<n;i++)
//	printf("%d ",a[i]);

    double Wr[MAX],Wi[MAX];
    Wr[1] = cos(-2. * M_PI / N);
    Wi[1] = sin(-2. * M_PI / N);
    Wr[0] = 1;
    Wi[0] = 0;
printf ("printing twiddle factors ....\n");
printf("1 %lf + %lfi\n", Wr[1], Wi[1]);
	
    for (int i = 2; i < N / 2; i++) {
        double angle = -2. * M_PI * i / N;
        Wr[i] = cos(angle);
        Wi[i] = sin(angle);
	printf("%d %lf , %lf\n", i , Wr[i],Wi[i]);
	//printf("%d %0b , %0b\n", i , Wr[i],Wi[i]);
    }
    printf("\n");

    int n = 1;
    int a = N / 2;

       // for(int i=0;i<N;i++)
//		printf("%lf + %lfi\n", fr[i], fi[i]);
//	printf("\n");

printf("log2%f\n",log2(16));
    for (int j = 0; j < log2(N); j++) {
        for (int i = 0; i < N; i++) {
            if (!(i & n)) {
                double  tempr = fr[i];
                double  tempi = fi[i];
                double Tempr,Tempi;
                Tempr = Wr[(i * a) % (n * a)] * fr[i + n] - Wi[(i * a) % (n * a)] * fi[i + n];
                Tempi = Wr[(i * a) % (n * a)] * fi[i + n] + Wi[(i * a) % (n * a)] * fr[i + n];
		printf("%lf + %lfi\n", tempr, tempi);
		printf("%lf + %lfi\n\n", Tempr, Tempi);
                fr[i] = tempr + Tempr;
                fi[i] = tempi + Tempi;
                fr[i + n] = tempr - Tempr;
                fi[i + n] = tempi - Tempi;
            }
        }
        n *= 2;
        a = a / 2;
    }


printf("...printing the FFT of the array specified\n");
for(int i=0 ;i<n;i++)
{
printf("%.5lf  + %.5lf j  \n",fr[i],fi[i]);
}

FILE *fp=fopen("nondsfftout.dat","w");
if(fp==NULL)
	printf("error file not found");
for (int i = 0; i < n; i++) {
    fprintf(fp, "%.5lf + %.5lfi\n", fr[i], fi[i]);
  }
fclose(fp);

}


void FFT(double fr[],double fi[],int n,double d){
transform(fr,fi,n);
for(int i=0;i<n;i++)
{
fr[i] *=d;
fi[i] *=d;
}
}


int main()
{
int a;
scanf("%d",&a);
//int n=131072;
int n=pow(2,a);
double d=1;

FILE *file;
double vecr[MAX],veci[MAX];

file =fopen("dummy.dat", "r" );
for (int i=0;i<n;i++)
{
veci[i]=0;
fscanf(file,"%lf",&vecr[i]);
}
fclose(file);

FFT(vecr,veci,n,d);

}

