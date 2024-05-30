/*
 * fix.h
 *
 * Fixed-point arithmetics.
 *
 * Using signed 32-bit integers: 12 bits for sign and integer part, 20 bits for fractional part
 *
 *  Created on: 12.10.2010
 *      Author: Alex Ignatov
 */

#ifndef __fix_H
#define __fix_H

#include <stdint.h>

typedef long long int int64_t;

#define fix_point 20				// bit length of fractional part. don't change, something will broke.
#define fix_base (1 << fix_point)
#define ONE (1 << fix_point)
#define PI 3294199					// just fixpoint PI constant

typedef int32_t fix_t;

#define _fix(x) ((int64_t)x << fix_point)

/* Simple conversion and arithmetic functions put in header file for further inlining */

static inline fix_t float2fix(float a) {
	return (fix_t)(a * fix_base);
}

static inline fix_t int2fix(int a) {
	return a << fix_point;
}

static inline float fix2float(fix_t a) {
	return a / (float)fix_base;
}

static inline fix_t fix_add(fix_t a, fix_t b) {
	return a + b;
}

static inline fix_t fix_sub(fix_t a, fix_t b) {
	return a - b;
}

/* Multiplication and division functions prefixes meaning:
 *
 * fix_fix_*		takes two fixed-point args, returns fixed-point
 * fix_int_*		takes fixed-point and integer args, returns fixed-point
 * int_int_*		takes two integer args, returns fixed-point
 */

static inline fix_t fix_fix_mul(fix_t a, fix_t b) {
	return (fix_t)(((int64_t)a * (int64_t)b) >> fix_point);
}

static inline fix_t fix_int_mul(fix_t a, int b) {
	return (fix_t)(a * b);
}

static inline fix_t fix_fix_div(fix_t a, fix_t b) {
	return (a < (0xFFFFFFFFU >> fix_point))
			? (fix_t)((a << fix_point) / b)
			: (fix_t)(((int64_t)a << fix_point) / b);
}

static inline fix_t int_int_div(int a, int b) {
	return (fix_t)(((int64_t)a << fix_point) / b);
}

static inline fix_t fix_int_div(fix_t a, int b) {
	return (fix_t)(a / b);
}

static inline fix_t fix_max(fix_t a, fix_t b) {
	return (a > b) ? a : b;
}

static inline fix_t fix_min(fix_t a, fix_t b) {
	return (a < b) ? a : b;
}

/* Shortcut macros for most frequently using functions */
#define fix_mul(a, b) fix_fix_mul(a, b)
#define fix_div(a, b) fix_fix_div(a, b)
#define fix_pow2(x) fix_mul(x, x)
#define fix_pow3(x) fix_mul(fix_pow2(x), x)
#define fix_pow4(x) fix_mul(fix_pow3(x), x)

#define fix_sign(a) (a & 0x80000000)
#define fix_abs(a) (fix_sign(a) ? -(a) : (a))

/* Trigonometry (via cordic) and other complex functions */
extern fix_t fix_atan2(fix_t y, fix_t x);
extern void fix_sincos(fix_t theta, fix_t *r_sin, fix_t *r_cos);
extern fix_t fix_sin(fix_t theta);
extern fix_t fix_cos(fix_t theta);
extern fix_t fix_tan(fix_t theta);
extern fix_t fix_sqrt(fix_t x);
extern fix_t fix_asin(fix_t x);
extern fix_t fix_acos(fix_t x);

#endif

