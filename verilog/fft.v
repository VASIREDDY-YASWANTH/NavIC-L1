module fft();
integer n=131072;
integer d=1;
integer signed arr [0:16];
integer i;


initial begin 
$readmemh("dummy.dat",arr);
for(i=0 ; i<16 ;i=i+1)
$display("array[%0d] = %b",i,arr[i]);
end

endmodule
