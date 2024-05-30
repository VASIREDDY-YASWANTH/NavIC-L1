module Testbench;

  reg signed [7:0] A[1:0][0:0];
  reg signed [7:0] B[1:0][0:0];
  reg signed [7:0] C[1:0][0:0];
  wire signed [7:0] sub[1:0][2:0];
  wire signed [7:0] add[1:0][2:0];
  wire signed [7:0] mid[1:0][2:0];
  wire signed [7:0] sub2[1:0][2:0];
  wire signed [7:0] med[1:0][2:0];
  wire signed [7:0] trans[2:0][1:0];
  wire signed [7:0] sc1[1:0][0:0];
  wire signed [7:0] sc2[1:0][0:0];
  wire signed [7:0] sc3[1:0][0:0];
  wire signed [7:0] sc4[1:0][0:0];
  wire signed [7:0] sc5[1:0][0:0];
  wire signed [7:0] sc6[1:0][0:0];
  wire signed [7:0] sc7[1:0][0:0];
  wire signed [7:0] sc8[1:0][0:0];
  wire signed [7:0] sc9[1:0][0:0];
  wire signed [7:0] N[1:0][1:0];
  wire signed [7:0] Q[1:0][1:0];
  wire signed [7:0] mul;
  wire signed [7:0] mul1;
  wire signed [7:0] Adet;
  wire signed [7:0] H[1:0][0:0];
  wire signed [7:0] med2[2:0][1:0];
  wire signed [7:0] sideAB;
  wire signed [7:0] sideBC;
  wire signed [7:0] sideCA;
  wire signed [7:0] x;
  wire signed [7:0] y;
  wire signed [7:0] z;
  wire signed [7:0] trans3[2:0][1:0];
  wire signed [7:0] xx;
  wire signed [7:0] yy;
  wire signed [7:0] zz;
  wire signed [7:0] ss;
  wire signed [7:0] dd;
  wire signed [7:0] ff;
  wire signed [7:0] l1;
  wire signed [7:0] l2;
  wire signed [7:0] l3;
  wire signed [15:0] fugito[1:0][0:0];
  wire signed [7:0] faggy[1:0][0:0];
  wire signed [7:0] foxy[1:0][0:0];
  wire signed [7:0] add2[1:0][0:0];
  wire signed [7:0] add3[1:0][0:0];
  wire signed [7:0] add34;
  wire signed [7:0] thank;
  wire signed [7:0] thank2;
  

  

 
  Matsub mysub (
    .A(A),
    .B(B),
    .C(C),
    .sc1(sc1),
    .sc2(sc2),
    .sc3(sc3),
    .sc4(sc4),
    .sc5(sc5),
    .sc6(sc6),
    .sub(sub),
    .add(add),
    .mid(mid),
    .sub2(sub2)
  );
  
  Matsub1 mysub1 (
    .A(A),
    .B(B),
    .C(C),
    .sub(sub),
    .med(med)
  );
  
  Mattran mattran(
    .A(A),
    .B(B),
    .C(C),
    .trans(trans)
  );

  Matscale matscale(
    .A(A),
    .B(B),
    .C(C),
    .l1(l1),
    .l2(l2),
    .l3(l3),
    .sc1(sc1),
    .sc2(sc2),
    .sc3(sc3),
    .sc4(sc4),
    .sc5(sc5),
    .sc6(sc6)
  );
       
  One one(
    .A(A),
    .B(B),
    .C(C),
    .sub(sub),
    .N(N)
  );   

  trans1 Trans1(
    .A(A),
    .B(B),
    .C(C),
    .N(N),
    .Q(Q)
  );   
  Matmul matmul(
     .A(A),
     .Q(Q),
     .B(B),
     .C(C),
     .mul(mul),
     .mul1(mul1)
  );
  Matdet matdet (
     .N(N),
     .mul(mul),
     .mul1(mul1),
     .Adet(Adet),
     .H(H)
 );
 
  huff uut(
       .med(med),
       .med2(med2)
 );
 
 side sides(
    .med(med),
    .med2(med2),
    .sideAB(sideAB),
    .sideBC(sideBC),
    .sideCA(sideCA)
 );
 
 hi norm(
    .sideAB(sideAB),
    .sideBC(sideBC),
    .sideCA(sideCA),
    .xx(xx),
    .yy(yy),
    .zz(zz),
    .x(x),
    .y(y),
    .z(z),
    .ss(ss),
    .dd(dd),
    .ff(ff)
);

  vatta pattakayibro (
         .A(A),
         .B(B),
         .C(C),
         .sub2(sub2),
         .trans3(trans3)
 );
 
 distance_mtr hibuddy (
      .trans3(trans3),
      .sub2(sub2),
      .xx(xx),
      .yy(yy),
      .zz(zz)
);

ignore_i capital( 
      .ss(ss),
      .dd(dd),
      .ff(ff),
      .l1(l1),
      .l2(l2),
      .l3(l3)
 );
 
Mat_add finesin (
     .add(add),
     .l1(l1),
     .l2(l2),
     .l3(l3),
     .fugito(fugito),
     .faggy(faggy),
     .foxy(foxy)
  );

 zum_m vroomm(
   .A(A),
   .B(B),
   .C(C),
   .ss(ss),
   .dd(dd),
   .ff(ff),
   .sc7(sc7),
   .sc8(sc8),
   .sc9(sc9)
 );
 
 fad_2 duffffer(
   .sc7(sc7),
   .sc8(sc8),
   .add2(add2)
 );
 
 fad_3 ffffer(
   .sc9(sc9),
   .add2(add2),
   .add3(add3)
 );
      
 nestedqt rofer (
     .ss(ss),
     .dd(dd),
     .ff(ff),
     .add34(add34)
 );
   
 tsal yfft (
    .add3(add3),
    .add34(add34),
    .thank(thank),
    .thank2(thank2)
);
   
   
    // Apply test vectors
  initial begin
    // Provide values for A, B, and C
    A[0][0] = 1; // Example value
    A[1][0] = -1; // Example value
    B[0][0] = -4; // Example value
    B[1][0] = 6; // Example value
    C[0][0] = -3; // Example value
    C[1][0] = -5; // Example value

    #1;
    $display("Input Values are: x1 = %d, y1 = %d, x2 = %d, y2 = %d, x3 = %d, y3 = %d", A[0][0], A[1][0], B[0][0], B[1][0], C[0][0], C[1][0]);

    #1;
    $display("MatSub value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 3; j = j + 1) begin
        $write("   %d", sub[i][j]);
      end
      $write("\n");
    end

    $display("Matadd value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 3; j = j + 1) begin
        $write("   %d", add[i][j]);
      end
      $write("\n");
    end

    $display("Matmid value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 3; j = j + 1) begin
        $write("   %d", mid[i][j]);
      end
      $write("\n");
    end

    $display("Matmed value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 3; j = j + 1) begin
        $write("   %d", med[i][j]);
      end
      $write("\n");
    end
  $display("Mattran value is:");
    #1;
    for (int i = 0; i < 3; i = i + 1) begin
      for (int j = 0; j < 2; j = j + 1) begin
        $write("   %d", trans[i][j]);
      end
      $write("\n");
    end
  $display("Matscale value is:");
#1;
for (int i = 0; i < 2; i = i + 1) begin
  for (int j = 0; j < 1; j = j + 1) begin
    $write("   sc1[%d][%d] = %d\n", i, j, sc1[i][j]);
    $write("   sc2[%d][%d] = %d\n", i, j, sc2[i][j]);
    $write("   sc3[%d][%d] = %d\n", i, j, sc3[i][j]);
    $write("   sc4[%d][%d] = %d\n", i, j, sc4[i][j]);
    $write("   sc5[%d][%d] = %d\n", i, j, sc5[i][j]);
    $write("   sc6[%d][%d] = %d\n", i, j, sc6[i][j]);
    $write("   sc7[%d][%d] = %d\n", i, j, sc7[i][j]);
    $write("   sc8[%d][%d] = %d\n", i, j, sc8[i][j]);
    $write("   sc9[%d][%d] = %d\n", i, j, sc9[i][j]);
    end
end

  $display("Matcomb value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 2; j = j + 1) begin
        $write("   %d", N[i][j]);
      end
      $write("\n");
    end
    $display("Mattran1 value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 2; j = j + 1) begin
        $write("   %d", Q[i][j]);
      end
      $write("\n");
    end
  $display("Matmul value is:");
    #1;
    $write("   %d", mul);
    $write("\n");
    $write("   %d", mul1);
    $write("\n");
    $write("   %d", Adet);
    $write("\n");
    $write(" %d ,%d", H[1][0] , H[0][0]);
    $write("\n");
    
    $display("Mattran2 value is:");
    #1;
    for (int i = 0; i < 3; i = i + 1) begin
      for (int j = 0; j < 2; j = j + 1) begin
        $write("   %d", med2[i][j]);
      end
      $write("\n");
    end
    $display("Matnorm value is:");
    #1;
    $write("   %d", sideAB);
    $write("\n");
    $write("   %d", sideBC);
    $write("\n");
    $write("   %d", sideCA);
    $write("\n"); 
    $write("   %d", x);
    $write("\n");
    $write("   %d", y);
    $write("\n");
    $write("   %d", z);
    $write("\n");
    $write("   %d", xx);
    $write("\n");
    $write("   %d", yy);
    $write("\n");
    $write("   %d", zz);
    $write("\n"); 
    $write("   %d", ss);
    $write("\n");
    $write("   %d",dd);
    $write("\n");
    $write("   %f", ff);
    $write("\n");
    $write("   %d", l1);
    $write("\n");
    $write("   %d", l2);
    $write("\n");
    $write("   %d", l3);
    $write("\n");
    $write(" hi  %d", add34);
    $write("\n");
    $write("   %d", thank);
    $write("\n");
    $write("   %d", thank2);
    $write("\n");
    
    
    $display("MatSub2 value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 3; j = j + 1) begin
        $write("   %d", sub2[i][j]);
      end
      $write("\n");
    end
    
    
    $display("Mattran3 value is:");
    #1;
    for (int i = 0; i < 3; i = i + 1) begin
      for (int j = 0; j < 2; j = j + 1) begin
        $write("   %d", trans3[i][j]);
      end
      $write("\n");
    end
    
   $display("Matscale3 value is:");
#1;
for (int i = 0; i < 2; i = i + 1) begin
  for (int j = 0; j < 1; j = j + 1) begin
    $write("   scale1[%d][%d] = %d\n", i, j, fugito[i][j]);
    $write("   scale2[%d][%d] = %d\n", i, j, faggy[i][j]);
    $write("   scale3[%d][%d] = %d\n", i, j, foxy[i][j]);
    end
end
     $display("Matadd3 value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 1; j = j + 1) begin
        $write("   %d", add2[i][j]);
        end
        $write("\n");
     end
  $display("Matadd4 value is:");
    #1;
    for (int i = 0; i < 2; i = i + 1) begin
      for (int j = 0; j < 1; j = j + 1) begin
        $write("   %d", add3[i][j]);
        end
        $write("\n");
     end
     
    
end    
endmodule





