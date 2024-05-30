module BCH;
  reg[51:0] encoded_symbols;
  reg feedback;
  reg[8:0] TOI_data;
  integer i;
  initial begin
    TOI_data = $urandom_range(400,1);
    //TOI_data = 9'd99;
    $display("TOI_data_int:%d",TOI_data);
    $display("TOI_data:%b",TOI_data);
    for (i=0; i<52; i=i+1) begin
      encoded_symbols[i] = TOI_data[8];
      feedback = TOI_data[8]^TOI_data[7]^TOI_data[6]^TOI_data[5]^TOI_data[4]^TOI_data[3]^TOI_data[1]^TOI_data[0];
      TOI_data[8:1] = TOI_data[7:0];
      TOI_data[0] = feedback;
    end
    $display("Encoded symbols:%b",encoded_symbols);
  end
endmodule
