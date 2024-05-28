
module fft();

function int Log2(int N);
static int k=N,i=0;
while(k) begin
k>>=1;
i++;
end

return i-1;
endfunction


function int reverse(int N,int n);
static int j,p=0;
for(j=1;j<=Log2(N);j++)
if(n&(1<< (Log2(N)-j)))
 p|=1<<(j-1);

return p;
endfunction


function int  transform(real fr[],real fi[],int N);
  
static int max=131072,n=1,a=N/2;
static real fr2[131072],fi2[131072],pi=4.0*$tan(1.0);
static real Wr[131072],Wi[131072],angle;
static real tempr,tempi,Tempr,Tempi;

for(int i=0;i<N;i++) 
    begin
      $display("a[%0d] = %.2f + i %.2f", i, fr[i],fi[i]);
    end
  $display("\n");

for (int i = 0; i < N; i++)
    begin
        fr2[i] = fr[reverse(N, i)];
        fi2[i] = fi[reverse(N, i)];
    end
    for (int j = 0; j < N; j++)
    begin
        fr[j] = fr2[j];
        fi[j] = fi2[j];
    end
 Wr[1] = $cos(-2.0 * pi / N);
    Wi[1] = $sin(-2.0 * pi / N);
    Wr[0] = 1;
    Wi[0] = 0;
     for (int i = 2; i < N / 2; i++) 
begin
        angle = -2.0 * pi * i / N;
        Wr[i] = $cos(angle);
        Wi[i] = $sin(angle);
end

 for (int j = 0; j < $clog2(N); j++) begin
        for (int i = 0; i < N; i++)
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

   for(int i=0;i<N;i++) 
    begin
      $display("a[%0d] = %.5f + i %.5f", i, fr[i],fi[i]);
    end


return N;
endfunction







   int n=16,d=1;
  real fr[16] = '{-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1};
  real fi[16] = '{0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
  initial begin
    for(int i=0;i<n;i++) 
    begin
   //   $display("a[%0d] = %.2f + i %.2f", i, fr[i],fi[i]);
    end
end
  
int x=Log2(n);
int z=transform(fr,fi,n);
initial begin 
$display("%0d",x);
end











endmodule
