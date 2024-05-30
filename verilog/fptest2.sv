

module testing2();

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




localparam N=16;
localparam SF = 2.0**-8.0;  // Q4.4 scaling factor is 2^-4
bit signed [15:0] fr[N],fi[N],fr2[N],fi2[N],tempr,tempi;  // (byte= 8bit 1 signed 3 integer 4 fractional)
bit signed [15:0] Wr[N],Wi[N] ,Tempr,Tempi;
int i,j,n=1,a=N/2;
int signed  ta,tb; // shortint =16 bit  using to store product of two 2 8bit values  

//reading data from files
initial begin
//$display("pi= %f",pi);

$readmemb("fr2.dat",fr);
$readmemb("fi2.dat",fi);
for (int i=0;i<N;i++)
$display("f[%0d] = %f +  %f j \n",i, $itor(fr[i]*SF), $itor(fi[i]*SF));
// seperating even and odd elements 
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

for (int i=0;i<N;i++)
$display("f[%0d] = %f +  %f j \n",i, $itor(fr[i]*SF), $itor(fi[i]*SF));

/*
//computing Twiddle factor
 	Wr[1] = $cos(-2.0 * pi / N);
 	Wi[1] = $sin(-2.0 * pi / N);
$display("\n%i",Wr[1]);
$display("\n%b",ab);
 	Wr[0] = 1;
 	Wi[0] = 0;
for (i = 2; i < N / 2; i++) 
	begin
        angle = -2.0 * pi * i / N;
        Wr[i] = $cos(angle);
        Wi[i] = $sin(angle);
	end
*/
// reading twiddle factor from file
$readmemb("wr2.dat",Wr);
$readmemb("wi2.dat",Wi);
$display("printing twiddle factors...");
for (int i=0;i<N/2;i++)
$display("f[%0d] = %f +  %f j \n",i, $itor(Wr[i]*SF), $itor(Wi[i]*SF));

//Cooley-Tukey FFT algorithm
for (int j = 0; j < Log2(N); j++) 
begin
  for ( int i = 0; i < N; i++)
  begin 
      if (!(i & n))
      begin
                tempr = fr[i];
                tempi = fi[i];
                ta=( Wr[(i * a) % (n * a)] * fr[i + n] );
		tb=( Wi[(i * a) % (n * a)] * fi[i + n] );
		Tempr=ta[15:8]-tb[15:8];     // considering only middle 8 bits of the product 
		ta= ( Wr[(i * a) % (n * a)] * fi[i + n] ); 
		tb= ( Wi[(i * a) % (n * a)] * fr[i + n] );
		Tempi=ta[15:8]+tb[15:8];	// considering only middle 8 bits of the product 
                fr[i] = tempr + Tempr;
                fi[i] = tempi + Tempi;
                fr[i + n] = tempr - Tempr;
                fi[i + n] = tempi - Tempi;
      end
   end    
        n *= 2;
        a = a/2;
end

//Printing  the FFT of input Array
$display("Printing the FFT of input Array...\n");
for(int i=0;i<N;i++) 
	begin
        $display("a[%0d] = %.5f + i %.5f", i, $itor(fr[i]*SF), $itor(fi[i]*SF));
	end

end 
endmodule 