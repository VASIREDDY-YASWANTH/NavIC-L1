module func_CORD();

// Declearing All variables used in code. 
localparam  ITERS=1024;
real theta_table[ITERS];
real pi=4.0*$atan(1.0);
real sin_x=0,cos_x=0,angle;
real K_n=1.0;
//real sigma,temp,theta,x,y,P2i,alpha;
int  j,n=ITERS;





// CORDIC algorithm 
function  real cord(real alpha);
real sigma,temp,theta,x,y,P2i;
 

 
	theta=0.0;x=1.0; y=0.0; P2i=1.0;
	  for( int i=0;i<n;i++)
	  begin
	  sigma=(theta<alpha)?1.0:-1.0;
	  theta+=sigma*theta_table[i];
	  //$display("theta=%.12f\n",theta);
	  temp=x;
	  x=x-sigma*y*P2i;
	  y=sigma*P2i*temp+y;
	  P2i/=2.0;
	  end

  //return 
 cos_x=x*K_n;
 sin_x=y*K_n;
endfunction



initial begin

// Generating Theta table.
   for (int i = 0; i < ITERS; i++)
   begin
   theta_table[i]=$atan2(1.0,$pow(2.0,i));
   end
// Computing K_n value.
   for( int  i=0;i<ITERS;i++) 
   begin
   K_n*=1.0/$sqrt(1.0+$pow(2.0,-2.0*i));
   end
//for (j = -90; j <= 90; j+=15)
  //begin
for (j = 2; j < 16; j++)
  begin
angle =-2*pi*j/32;
//angle= j*pi/180;
cord(angle);
  $display("degree     sin_x       diff");
  $display("%0d %.8f (%.16f) \n",j,sin_x,  (sin_x-$sin(angle) ) );
  end

end
endmodule
