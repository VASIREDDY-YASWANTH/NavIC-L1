
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

initial begin
 for(int i=0;i<N;i++) 
    begin
   //   $display("a[%0d] = %.2f + i %.2f", i, fr[i],fi[i]);
    end
end


initial begin

for(i=0;i<N;i++) 
    begin
      $display("a[%0d] = %.2f + i %.2f", i, fr[i],fi[i]);
    end
  $display("\n");
 


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
 for( i=0;i<N;i++) 
    begin
      $display("a[%0d] = %.5f + i %.5f", i, fr[i],fi[i]);
    end
$display("\n");

 Wr[1] = $cos(-2.0 * pi / N);
    Wi[1] = $sin(-2.0 * pi / N);
$display(" wr1 =%.5f ,wi1 =%.5f\n",  Wr[1],Wi[1]);
    Wr[0] = 1;
    Wi[0] = 0;
////////////////
     for (i = 2; i < N / 2; i++) 
begin
        angle = -2.0 * pi * i / N;
        Wr[i] = $cos(angle);
        Wi[i] = $sin(angle);
$display(" wr =%.5f ,wi =%.5f\n",  Wr[i],Wi[i]);
end
////////////
    for (j = 0; j < N; j++)
    begin
        ftr[j] = fr[j];
        fti[j] = fi[j];
    end

for (int j = 0; j < Log2(N); j++) 
begin
  for ( int i = 0; i < N; i++)
  begin 
      if (!(i & n))
      begin
                tempr = ftr[i];
                tempi = fti[i];
                Tempr = Wr[(i * a) % (n * a)] * ftr[i + n] - Wi[(i * a) % (n * a)] * fti[i + n];
                Tempi = Wr[(i * a) % (n * a)] * fti[i + n] + Wi[(i * a) % (n * a)] * ftr[i + n];
	$display(" tempr =%.5f ,tempi =%.5f\n", tempr,tempi);
        $display(" Tempr =%.5f ,Tempi =%.5f\n",  Tempr,Tempi);
                ftr[i] = tempr + Tempr;
                fti[i] = tempi + Tempi;
                ftr[i + n] = tempr - Tempr;
                fti[i + n] = tempi - Tempi;
       end
  end    
        n *= 2;
        a = a / 2;
 end


   for(int i=0;i<N;i++) 
    begin
      $display("a[%0d] = %.5f + i %.5f", i, ftr[i],fti[i]);
    end
end



endmodule
