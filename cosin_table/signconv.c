#include <stdio.h>

#define INT_BITS 16  // Number of integer bits
#define FRAC_BITS 16 // Number of fractional bits

// Function to convert integer to fixed-point
int intToFixedPoint(int num) {
    return num << FRAC_BITS;
}

int main() {
    int integerNumber = 42;
    int fixedPointNumber = intToFixedPoint(integerNumber);

    printf("Integer: %d\nFixed-Point: %d\n", integerNumber, fixedPointNumber);

    return 0;
}
