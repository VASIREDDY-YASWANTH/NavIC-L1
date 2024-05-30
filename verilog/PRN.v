module PRN_generator;
  reg[0:54] R0;
  reg[0:54] R1;
  reg[0:4] C;
  reg[0:10229] PRN;
  integer i;
  integer j;
  integer p;
  reg R0_fb;
  reg sigma_2A;
  reg sigma_2B;
  reg sigma_2C;
  reg sigma_2;
  reg R1A;
  reg R1B;
  reg R1_fb;
  reg C_fb;
  initial begin
    R0 = 55'o0227743641272102303;
    R1 = 55'o1667217344450257245;
    C = 5'b01000;
    for(i=0; i<10230;i=i+1) begin
      R0_fb = R0[50]^R0[45]^R0[40]^R0[20]^R0[10]^R0[5]^R0[0];
      sigma_2A = (R0[50]^R0[45]^R0[40])&(R0[20]^R0[10]^R0[5]^R0[0]);
      sigma_2B = ((R0[50]^R0[45])&R0[40])^((R0[20]^R0[10])&(R0[5]^R0[0]));
      sigma_2C = (R0[50]&R0[45])^(R0[20]&R0[10])^(R0[5]&R0[0]);
      sigma_2 = sigma_2A^sigma_2B^sigma_2C;
      R1A = sigma_2 ^ (R0[40]^R0[35]^R0[30]^R0[25]^R0[15]^R0[0]);
      R1B = R1[50]^ R1[45]^R1[40]^ R1[20]^ R1[10]^ R1[5]^ R1[0];
      R1_fb = R1A ^ R1B;
      PRN[i] = R1[0]^C[0];
      R0[0:53] =  R0[1:54]; 
      R0[54] = R0_fb;
      R1[0:53] = R1[1:54];
      R1[54] = R1_fb;
      C_fb = C[0];
      C[0:3] = C[1:4];
      C[4] = C_fb;
   end
    $display("First 24: %o",PRN[0:23]);
    $display("Last 24: %o",PRN[10206:10229]);
  end
endmodule
