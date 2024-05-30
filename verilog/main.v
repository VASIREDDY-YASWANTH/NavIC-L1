module Matsub (
  input signed [7:0] A[1:0][0:0],
  input signed [7:0] B[1:0][0:0],
  input signed [7:0] C[1:0][0:0],
  input signed [7:0] sc1[1:0][0:0],
  input signed [7:0] sc2[1:0][0:0],
  input signed [7:0] sc3[1:0][0:0],
  input signed [7:0] sc4[1:0][0:0],
  input signed [7:0] sc5[1:0][0:0],
  input signed [7:0] sc6[1:0][0:0],
  output reg signed [7:0] sub[1:0][2:0],
  output reg signed [7:0] add[1:0][2:0],
  output reg signed [7:0] mid[1:0][2:0],
  output reg signed [7:0] sub2[1:0][2:0]
 );

  integer i,j;
  
  always @* begin
    for (i = 0; i < 2; i = i+1) begin
        sub[i][0] = A[i][0] - B[i][0];
        sub[i][1] = B[i][0] - C[i][0];
        sub[i][2] = C[i][0] - A[i][0];
        add[i][0] = sc1[i][0] + sc2[i][0];
        add[i][1] = sc3[i][0] + sc4[i][0];
        add[i][2] = sc5[i][0] + sc6[i][0];
        mid[i][0] = (A[i][0] + B[i][0])/2;
        mid[i][1] = (B[i][0] + C[i][0])/2;
        mid[i][2] = (C[i][0] + A[i][0])/2;
        sub2[i][0] = B[i][0] - A[i][0];
        sub2[i][1] = C[i][0] - A[i][0];
        sub2[i][2] = C[i][0] - B[i][0];
        
    end
  end
endmodule
module Matsub1(
  input signed [7:0] A[1:0][0:0],
  input signed [7:0] B[1:0][0:0],
  input signed [7:0] C[1:0][0:0],
  input signed [7:0] sub[1:0][2:0], // Input sub array from Matsub
  output reg signed [7:0] med[1:0][2:0]
);
  integer i,j;
  
  always @* begin
    for (int i = 0; i < 3; i = i + 1) begin // Corrected loop limit
        med[i][0] = A[i][0] - sub[i][0];
        med[i][1] = B[i][0] - sub[i][1];
        med[i][2] = C[i][0] - sub[i][2];
    end
  end
endmodule

module Mattran(
  input signed [7:0] A[1:0][0:0],
  input signed [7:0] B[1:0][0:0],
  input signed [7:0] C[1:0][0:0],
  output reg signed [7:0] trans[2:0][1:0]
);
integer i, j;
always @* begin
 for(i=0;i<3;i=i+1) begin
  for(j=0;j<2;j=j+1) begin
    trans[0][j] = A[j][0];
    trans[1][j] = B[j][0];
    trans[2][j] = C[j][0]; 
  end
 end
end
endmodule

  
module Matscale(
  input signed [7:0] A[1:0][0:0],
  input signed [7:0] B[1:0][0:0],
  input signed [7:0] C[1:0][0:0],
  input signed [7:0] l1,
  input signed [7:0] l2, // Resulting 1x1 matrix
  input signed [7:0] l3,
  output reg signed [7:0] sc1[1:0][0:0],
  output reg signed [7:0] sc2[1:0][0:0],
  output reg signed [7:0] sc3[1:0][0:0],
  output reg signed [7:0] sc4[1:0][0:0],
  output reg signed [7:0] sc5[1:0][0:0],
  output reg signed [7:0] sc6[1:0][0:0]
  );
 integer i,j;
 always @* begin
   for(i=0;i<2;i++) begin
  sc1[i][0]=l1*C[i][0];
  sc2[i][0]=l2*B[i][0];
  sc3[i][0]=l2*A[i][0];
  sc4[i][0]=l3*C[i][0];
  sc5[i][0]=l3*B[i][0];
  sc6[i][0]=l1*A[i][0];
  end
 end 
endmodule    
  
module One(
    input signed [7:0] A[1:0][0:0],
    input signed [7:0] B[1:0][0:0],
    input signed [7:0] C[1:0][0:0],
    input signed [7:0] sub[1:0][2:0], // Input sub array from Matsub
    output reg signed [7:0] N[1:0][1:0] // Output transposed matrix
);
integer i;
always @* begin
    for (i = 0; i < 2; i = i + 1) begin
        N[i][0] = 0+sub[i][0]; // Extract sub[i][0]
        N[i][1] = 0+sub[i][2]; // Extract sub[i][2]
    end
end
endmodule

module trans1(
       input signed [7:0] A[1:0][0:0],
       input signed [7:0] B[1:0][0:0],
       input signed [7:0] C[1:0][0:0],
       input signed [7:0] N[1:0][1:0], // Input sub array from Matsub
       output reg signed [7:0] Q[1:0][1:0]
);
 integer i, j;
always @* begin
 for(i=0;i<3;i=i+1) begin
  for(j=0;j<2;j=j+1) begin
    Q[0][j] = N[j][0];
    Q[1][j] = N[j][1];
    end
 end
end
endmodule

module Matmul(
    input signed [7:0] A[1:0][0:0],
    input signed [7:0] Q[1:0][1:0],
    input signed [7:0] B[1:0][0:0],
    input signed [7:0] C[1:0][0:0],
    output reg signed [7:0] mul, // Resulting 1x1 matrix
    output reg signed [7:0] mul1
);

reg signed [15:0] temp; // Temporary variable for storing intermediate result
reg signed [15:0] temp1;

always @* begin
    temp=0;
    temp1 = 0; // Initialize temporary variable to zero
    for (int i = 0; i < 2; i = i + 1) begin
        temp = temp + (Q[1][i] * B[i][0]); // Perform matrix multiplication
        temp1 = temp1+(Q[0][i] * C[i][0]);
    end
    mul = temp;
    mul1=temp1; // Store the result in the output matrix
end
endmodule

module Matdet(
      input signed [7:0] N[1:0][1:0],
      input signed [7:0] mul, // Resulting 1x1 matrix
      input signed [7:0] mul1,
      output reg signed [7:0] Adet,
      output reg signed [7:0] H[1:0][0:0]
   );
 reg signed [15:0] temp;     
      always @*begin
      Adet = N[0][0]*N[1][1]-N[0][1]*N[1][0];
      temp=Adet;
      H[1][0]=(-1)*(mul*N[1][0]-N[0][0]*mul1)/temp;
      H[0][0]=(mul*N[1][1]-N[0][1]*mul1)/temp;
      end 
 endmodule  


module huff(
   input signed [7:0] med[1:0][2:0],
   output reg signed [7:0] med2[2:0][1:0]
  ); 
   integer i, j;
  always @* begin
  for(i=0;i<3;i=i+1) begin
  for(j=0;j<2;j=j+1) begin
    med2[0][j] = 0+med[j][0];
    med2[1][j] = 0+med[j][1];
    med2[2][j] = 0+med[j][2]; 
  end
 end
end
endmodule

module side(
      input signed [7:0] med[1:0][2:0],
      input signed [7:0] med2[2:0][1:0],
   output reg signed [7:0] sideAB,
   output reg signed [7:0] sideBC, // Resulting 1x1 matrix
    output reg signed [7:0] sideCA
);

reg signed [15:0] temp; // Temporary variable for storing intermediate result
reg signed [15:0] temp1;
reg signed [15:0] temp2;

always @* begin
    temp=0;
    temp1 = 0; // Initialize temporary variable to zero
    temp2=0;
    for (int i = 0; i < 2; i = i + 1) begin
        temp = temp + (med2[0][i] * med[i][0]);
        temp1 = temp1+(med2[1][i] * med[i][1]);
        temp2 = temp2+(med2[2][i] * med[i][2]);
    
    end
    sideAB = temp;
    sideBC=  temp1;
    sideCA=  temp2; // Store the result in the output matrix
 end
endmodule

module hi(
      input signed [7:0] sideAB,
      input signed [7:0] sideBC, // Resulting 1x1 matrix
      input signed [7:0] sideCA,
      input signed [7:0] xx,
      input signed [7:0] yy, // Resulting 1x1 matrix
      input signed [7:0] zz,
      output reg signed [7:0] x,
      output reg signed [7:0] y, // Resulting 1x1 matrix
      output reg signed [7:0] z,
      output reg signed [7:0] ss,
      output reg signed [7:0] dd, // Resulting 1x1 matrix
      output reg signed [7:0] ff
   );
      always @*begin
      
      x= $sqrt(sideAB);
      y= $sqrt(sideBC);
      z= $sqrt(sideCA);
      ss= $sqrt(xx);
      dd= $sqrt(yy);
      ff= $sqrt(zz);
 end
endmodule

module vatta(
   input signed [7:0] A[1:0][0:0],
   input signed [7:0] B[1:0][0:0],
   input signed [7:0] C[1:0][0:0],
   input signed [7:0] sub2[1:0][2:0],
   output reg signed [7:0] trans3[2:0][1:0]
  ); 
   integer i, j;
  always @* begin
  for(i=0;i<3;i=i+1) begin
  for(j=0;j<2;j=j+1) begin
    trans3[0][j] = sub2[j][0];
    trans3[1][j] = sub2[j][1];
    trans3[2][j] = sub2[j][2]; 
  end
 end
end
endmodule

module distance_mtr(
      input signed [7:0] sub2[1:0][2:0],
      input signed [7:0] trans3[2:0][1:0],
   output reg signed [7:0] xx,
   output reg signed [7:0] yy, // Resulting 1x1 matrix
    output reg signed [7:0] zz
);

reg signed [15:0] temp; // Temporary variable for storing intermediate result
reg signed [15:0] temp1;
reg signed [15:0] temp2;

always @* begin
    temp=0;
    temp1 = 0; // Initialize temporary variable to zero
    temp2=0;
    for (int i = 0; i < 2; i = i + 1) begin
        temp = temp + (trans3[0][i] * sub2[i][0]);
        temp1 = temp1+(trans3[1][i] * sub2[i][1]);
        temp2 = temp2+(trans3[2][i] * sub2[i][2]);
    
    end
    xx = temp;
    yy=  temp1;
    zz=  temp2;
     // Store the result in the output matrix
 end
endmodule

module ignore_i (
     input signed [7:0] ss,
      input signed [7:0] dd, // Resulting 1x1 matrix
      input signed [7:0] ff,
      output reg signed [7:0] l1,
      output reg signed [7:0] l2, // Resulting 1x1 matrix
      output reg signed [7:0] l3
 );
  
  always @* begin
  
  l1=(ff+ss-dd)/2;
  l2=(ff+dd-ss)/2;
  l3=(ss+dd-ff)/2;
  
  end
endmodule 

module Mat_add(
  input signed [7:0] add[1:0][2:0],
  input signed [7:0] l1,
  input signed [7:0] l2, // Resulting 1x1 matrix
  input signed [7:0] l3,
  output reg signed [15:0] fugito[1:0][0:0],
  output reg signed [7:0] faggy[1:0][0:0],
  output reg signed [7:0] foxy[1:0][0:0]
);
 always @* begin

  for(int i=0;i<2;i++) begin
  fugito[i][0]=1/(l1+l2) * add[i][0];
  faggy[i][0]=1/(l2+l3) * add[i][1];
  foxy[i][0]=1/(l3+l1) * add[i][2];

end
 end 
endmodule 

module zum_m (
     input signed [7:0] ss,
      input signed [7:0] dd, // Resulting 1x1 matrix
      input signed [7:0] ff,
      input signed [7:0] A[1:0][0:0],
      input signed [7:0] B[1:0][0:0],
      input signed [7:0] C[1:0][0:0],
      output reg signed [7:0] sc7[1:0][0:0],
      output reg signed [7:0] sc8[1:0][0:0],
      output reg signed [7:0] sc9[1:0][0:0]
 );
 
   always @* begin
   for(int i=0;i<2;i++) begin
  sc7[i][0]=ff*A[i][0];
  sc8[i][0]=dd*B[i][0];
  sc9[i][0]=ss*C[i][0];
 end
end
endmodule

module fad_2 (
  input signed [7:0] sc7[1:0][0:0],
  input signed [7:0] sc8[1:0][0:0],
 output reg signed [7:0] add2[1:0][0:0]
  
  );
  
  integer i;
  
  always @* begin
    for (i = 0; i < 2; i = i+1) begin
        add2[i][0] = sc7[i][0] + sc8[i][0];
        end
 end
endmodule      

module fad_3 (
  input signed [7:0] sc9[1:0][0:0],
  input signed [7:0] add2[1:0][0:0],
 output reg signed [7:0] add3[1:0][0:0]
  
  );
  
  integer i;
  
  always @* begin
    for (i = 0; i < 2; i = i+1) begin
        add3[i][0] = add2[i][0] + sc9[i][0];
        end
 end
endmodule  

module nestedqt (
      input signed [7:0] ss,
      input signed [7:0] dd, // Resulting 1x1 matrix
      input signed [7:0] ff,
      output reg signed [7:0] add34
 );
 
 always @* begin
    
    add34 = ff+dd+ss;
  end
 endmodule
 
module tsal (
       input signed [7:0] add3[1:0][0:0],
       input signed [7:0] add34,
       output reg signed [7:0] thank,
       output reg signed [7:0] thank2
       );

always @* begin
  
    thank=add3[0][0]/add34;
    thank2=add3[1][0]/add34;
  end
endmodule
    
     
       
      
 



    
    



   
   
     
   





