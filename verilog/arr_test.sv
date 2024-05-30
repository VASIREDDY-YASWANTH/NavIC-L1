module arr;
  bit signed [15:0] ar[16] = '{16'b0_000000000000000, 16'b0_000110010001011, 16'b0_001100011111001, 16'b0_010010100101000, 16'b0_011000011111100, 16'b0_011110001010110, 16'b0_100011100011101, 16'b0_101000100110100, 16'b0_101101010000010, 16'b0_110001011110001, 16'b0_110101001101110, 16'b0_111000011100010, 16'b0_111011001000010, 16'b0_111101001111011, 16'b0_111110110001010, 16'b0_111111101100010
};
  localparam SF = 2.0**-15.0;
  int n = 1;
  function bit [15:0] sin_tab;
    input [3:0] ang;
    begin
      sin_tab = ar[ang];
    end
  endfunction
  initial begin
    bit [3:0] val;
    bit signed [15:0]result;
    val =11;
    result = sin_tab(val);
    $display("array[%0d] = %0f", val, result*SF);
  end
endmodule