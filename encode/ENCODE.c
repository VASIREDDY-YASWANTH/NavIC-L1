#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include"NavicL1.h"
#include <inttypes.h> // Include this header for uint32_t type
int main()
{

int r=1 ,c=250;
// load buffer array
    char ***buff = loadStrArr("buff2.dat", r, c, 1); 
printStrArr(buff,r,c);

}
