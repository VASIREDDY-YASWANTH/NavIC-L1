/*
 * fix.c
 *
 * Fixed-point arithmetics.
 *
 *  Created on: 12.10.2010
 *      Author: Alex Ignatov
 */

#include "fix.h"

#include "cos_table.h"

#define CORDIC_COUNT 16
const fix_t cordic[] = { 823549, 486169, 256878, 130395, 65450, 32757, 16382, 8191, 4095, 2047, 1023, 511, 255, 127, 63, 31 };
#define K1 (int32_t)(0.6072529350088812561694 * fix_base)

#define HALF_PI PI/2

/* atan2 using CORDIC method */
fix_t fix_atan2(fix_t y, fix_t x) {
	int32_t theta, i, t;

	theta = 0;
	if (x<0) {
		if (y >= 0) {
			theta = HALF_PI;
			t = x;
			x = y;
			y = -t;
		}
		else {
			theta = -HALF_PI;
			t = x;
			x = -y;
			y = t;
		}
	}
	else
		theta = 0;

	int32_t last_x;
	for (i = 0; i < CORDIC_COUNT; i++) {
		last_x = x;

		if (y < 0) { // sign=1
			x -= y >> i;
			y += last_x >> i;
			theta -= cordic[i];
		}
		else {
			x += y >> i;
			y -= last_x >> i;
			theta += cordic[i];
		}
   }
   return theta;
}

void fix_sincos(fix_t theta, fix_t *r_sin, fix_t *r_cos) {
	if ((theta >= COS_MIN) && (theta < COS_MAX)) {
		// angle small enough to use linear approximation for sin and table lookup for cos
		*r_sin = theta;
		*r_cos = cos_table[(theta - COS_MIN) >> COS_SHIFT];
		return;
	}

	int32_t i, d, tx, ty, tz;
	int32_t x = K1, y = 0, z;
	uint8_t flip = 0; // flip sign

	// rotate argument to I quadrant, remembering sign changes
	if (theta < -HALF_PI || theta > HALF_PI) {
		if (theta < 0)
			z = theta + PI;
		else
			z = theta - PI;
		flip = 1;
	}
	else {
		z = theta;
		flip = 0;
	}

	// actual CORDIC algorithm, see wiki
	for (i = 0; i < CORDIC_COUNT; i++) {
	    d = z>>31;
    	tx = x - (((y>>i) ^ d) - d);
		ty = y + (((x>>i) ^ d) - d);
		tz = z - ((cordic[i] ^ d) - d);
    	x = tx; y = ty; z = tz;
	}
	if (flip) {
		x = -x;
		y = -y;
	}
	*r_cos = x;
	*r_sin = y;
}

/* fix_sin, fix_cos, fix_tan - wrappers for fix_sincos */
fix_t fix_sin(fix_t theta) {
	fix_t r_sin, r_cos;
	fix_sincos(theta, &r_sin, &r_cos);
	return r_sin;
}

fix_t fix_cos(fix_t theta) {
	fix_t r_sin, r_cos;
	fix_sincos(theta, &r_sin, &r_cos);
	return r_cos;
}

fix_t fix_tan(fix_t theta) {
	fix_t r_sin, r_cos;
	fix_sincos(theta, &r_sin, &r_cos);
	return fix_div(r_sin, r_cos);
}

/* Square root using cordic-like algorithm
 * http://www.convict.lu/Jeunes/Math/square_root_CORDIC.htm
 */
fix_t fix_sqrt(fix_t x) {
	fix_t base, y;
	int i;
	if (x == 0)
		return 0;
	base = 33554432;
	y = 0;
	for (i = 0; i < 26; i++) {
		y += base;
		if  (fix_mul(y, y) > x)
			y -= base;
		base >>= 1;
	}
	return y;
}

fix_t fix_asin(fix_t x) {
	// asin(x) = atan(x/sqrt(1-x^2))
	return fix_atan2(x, fix_sqrt(ONE - fix_pow2(x)));
}

fix_t fix_acos(fix_t x) {
	// acos(x) = atan(sqrt(1-x^2)/x)
	return fix_atan2(fix_sqrt(ONE - fix_pow2(x)), x);
}

