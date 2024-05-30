
module fft();

// function for log base2
function int Log2(int N);
static int k,i;
k=N;i=0;
while(k) begin
k>>=1;
i++;
end
return i-1;
endfunction

//function to seperate even and odd elements
function int reverse(int N,int n);
static int j,p;
p=0;
for(j=1;j<=Log2(N);j++)
if(n&(1<< (Log2(N)-j)))
 p|=1<<(j-1);
return p;
endfunction

//declearing All variables used in code 
int N=32,d=1;
real fr[32] = '{-1,-1,-1,-1, -1,-1,1,1, 1,1,1,1, 1,1,1,1, -1,-1,-1,-1, -1,1,1,1, 1,1,-1,-1, -1,-1,-1,-1};
real fi[32] = '{0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
int max=131072,n=1,a=N/2,i=0,j=0;
real fr2[131072],fi2[131072],ftr[131072],fti[131072],pi=4.0*$atan(1.0);
real Wr[131072],Wi[131072],angle;
real tempr,tempi,Tempr,Tempi;

// FFT Algorithm
initial begin

for ( i = 0; i < N; i++)
	begin
        fr2[i] = fr[reverse(N, i)];
        fi2[i] = fi[reverse(N, i)];
	end
for (j = 0; j < N; j++)
	begin
        fr[j] = fr2[j];
        fi[j] = fi2[j];
	end

//computing Twiddle factor
 	Wr[1] = $cos(-2.0 * pi / N);
 	Wi[1] = $sin(-2.0 * pi / N);
 	Wr[0] = 1;
 	Wi[0] = 0;
for (i = 2; i < N / 2; i++) 
	begin
        angle = -2.0 * pi * i / N;
        Wr[i] = $cos(angle);
        Wi[i] = $sin(angle);
	end
//Cooley-Tukey FFT algorithm
for (int j = 0; j < Log2(N); j++) 
begin
  for ( int i = 0; i < N; i++)
  begin 
      if (!(i & n))
      begin
                tempr = fr[i];
                tempi = fi[i];
                Tempr = Wr[(i * a) % (n * a)] * fr[i + n] - Wi[(i * a) % (n * a)] * fi[i + n];
                Tempi = Wr[(i * a) % (n * a)] * fi[i + n] + Wi[(i * a) % (n * a)] * fr[i + n];
                fr[i] = tempr + Tempr;
                fi[i] = tempi + Tempi;
                fr[i + n] = tempr - Tempr;
                fi[i + n] = tempi - Tempi;
      end
   end    
        n *= 2;
        a = a / 2;
end

//Printing  the FFT of input Array
$display("Printing the FFT of input Array...\n");
for(int i=0;i<N;i++) 
	begin
        $display("a[%0d] = %.5f + i %.5f", i, fr[i],fi[i]);
	end

end

endmodule