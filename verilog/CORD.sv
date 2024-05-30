module CORDIC();

// Declearing All variables used in code. 
localparam  ITERS=32;
real theta_table[ITERS];
real pi=4.0*$atan(1.0);
real sin_x=0;
real K_n=1.0;
real sigma,temp,theta,x,y,P2i,alpha;
int  i,j,n;

// CORDIC algorithm 
initial begin
   // Generating Theta table.
   for (i = 0; i < ITERS; i++)
   begin
   theta_table[i]=$atan2(1.0,$pow(2.0,i));
   end
   // Computing K_n value.
   for( i=0;i<ITERS;i++) 
   begin
   K_n*=1.0/$sqrt(1.0+$pow(2.0,-2.0*i));
   end

  for (j = -90; j <= 90; j+=15)
  begin
	theta=0.0;x=1.0; y=0.0; P2i=1.0;
	alpha=j*(pi/180);
	n=ITERS;
	  for( i=0;i<n;i++)
	  begin
	  sigma=(theta<alpha)?1.0:-1.0;
	  theta+=sigma*theta_table[i];
	  //$display("theta=%.12f\n",theta);
	  temp=x;
	  x=x-sigma*y*P2i;
	  y=sigma*P2i*temp+y;
	  P2i/=2.0;
	  end
  sin_x= y*K_n;
  $display("degree     sin_x       diff");
  $display("%0d %.8f (%.8f) \n",j,sin_x,  (sin_x-$sin(j*(pi/180)) ) );
  end

end
endmodule