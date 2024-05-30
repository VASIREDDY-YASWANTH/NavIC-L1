module fixedtest();
  bit signed [7:0] a, b, c;
  bit signed [15:0] ab;  // large enough for product
  int d,fs=18548387,length=10230,basis=1023000;
  bit signed [63:0] count,f,g;  

  localparam SF = 2.0**-4.0;
  localparam SF2 = 2.0**-32.0; // Q4.4 scaling factor is 2^-4
  localparam SF3 = 2.0**-8.0;
  
  initial begin
//     $display("Fixed Point Examples from projectf.io.");

//     a = 8'b0011_1010;  // 3.6250
//     b = 8'b0100_0001;  // 4.0625
//     c = a + b;         // 0111.1011 = 7.6875
//     $display("%f + %f = %f", a, $itor(b*SF), $itor(c*SF));

//     a = 8'b0011_1010;  // 3.6250
//     b = 8'b1110_1000;  // -1.5000
//     c = a + b;         // 0010.0010 = 2.1250
//     $display("%f + %f = %f", $itor(a*SF), $itor(b*SF), $itor(c*SF));

//     a = 8'b0011_0100;  // 3.2500
//     b = 8'b0010_0001;  // 2.0625
//     ab = a * b;        // 00000110.10110100 = 6.703125
//     c = ab[11:4];      // take middle 8 bits: 0110.1011 = 6.6875
//     $display("%f * %f = %f", $itor(a*SF), $itor(b*SF), $itor(c*SF));

    a = 8'b1010_1000;  // 7.5000
    b = 8'b0111_0110;  // 0.5000
    ab = a * b;        // 00000011.11000000 = 3.7500
    d = a * b;
    f = fs;
    g = (f<<32);
    count = (g/(basis/length))*SF2;
    $display("%b",ab);
    $display("%b",d);
    $display("%b",f);
    $display("%b",g);
    //c = a*(2.0**-5.0);      // take middle 8 bits: 0011.1100 = 3.7500
    //c = a;
    $display("%f * %f = %g", $itor(a*SF), $itor(b*SF), $itor(ab*SF3));
    $display("Result: %g",$itor(d*SF3));
    $display("%f",count);
  end
endmodule
