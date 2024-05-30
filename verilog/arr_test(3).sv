module arr;
  `define pi 16'b011_0010010000111
  bit signed [15:0] array[16] = '{16'b0_000000000000000, 16'b0_000110010001011, 16'b0_001100011111001, 16'b0_010010100101000, 16'b0_011000011111100, 16'b0_011110001010110, 16'b0_100011100011101, 16'b0_101000100110100, 16'b0_101101010000010, 16'b0_110001011110001, 16'b0_110101001101110, 16'b0_111000011100010, 16'b0_111011001000010, 16'b0_111101001111011, 16'b0_111110110001010, 16'b0_111111101100010
                              };
  localparam SF = 2.0**-15.0;
  localparam SF1 = 2.0**-13.0;
  int n = 1;
  function bit signed [15:0] sin_tab;
    input bit [15:0] angle; 
    int index;
    begin
index = int'(angle/(`pi/(2*16)));
$display("Index = %0d",index);
if (index<16)begin
      $display("Index = %0d",index);
      sin_tab = array[index];
end
 else if (32>index)begin
        index = 31-index;
      $display("Index = %0d",index);
      sin_tab = array[index];
      end
else if (48>index)begin
        index = index-32;
 $display("Index = %0d",index);
      sin_tab = -array[index];
      end
else if (64>index)begin
        index = 63-index;
 $display("Index = %0d",index);
      sin_tab = -array[index];
      end

    end
  endfunction
  initial begin
    bit [15:0] val;
    bit signed [15:0]result;
    val =5*(`pi/3);
    result = sin_tab(val);
    $display("sin(%0f) = %0f", val*SF1, result*SF);
  end
endmodule