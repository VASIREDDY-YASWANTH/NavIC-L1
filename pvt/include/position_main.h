/* header file for posiition_main file */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>


typedef int                int_32;
typedef unsigned int       uint_32;
typedef unsigned char      uint_8;
typedef signed char        int_8;
typedef unsigned long long uint_64;

static struct Space_vehicle{
    int ids;
    double mean_motion_at_ephemeris; // M0
    double sqrt_A;
    double Cus;
    double Cuc;
    double Cis;
    double Cic;
    double Crc;
    double Crs;
    double delta_n;
    double e;
    double w;
    double omega_0;
    double idot;
    double omega_dot;
    double inclination_at_ephemeris;
    double group_delay;
    double reference_time;
    double correction_f2;
    double correction_f1;
    double correction_f0;
    double pos_t;
    double time_correction;
    double time_of_ephemeris;
    uint_32 lock_code_nco;
    double time_raw;
    double nav_subframe_week;
    double lock_ms_frame;
    uint_32 subframe[6][11];
    uint_8 subframe_number;
    uint_8 synced;
}space_vehicles[] = {{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32}};

struct Location {
    double x;
    double y;
    double z;
    double time;
};

void nav_save_frame(struct Space_vehicle sv);
void sv_calc_corrected_time(int id);
double  orbit_ecc_anom(struct Space_vehicle *sv, double t);
int sv_calc_location(int id, struct Location *l);
int orbit_calc_position(struct Space_vehicle *sv, struct Location *l);
int A_solve(int chans, struct Location *sv_l, struct Location *p);
void LatLonAlt( double x_n, double y_n, double z_n, double lat, double lon, double alt);
