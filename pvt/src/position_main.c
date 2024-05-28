 /* main function that calculates receiver location */
#include "position_main.h"

#define MAX_ITER 20
#define WGS84_A     (6378137.0)
#define WGS84_F_INV (298.257223563)
#define WGS84_B     (6356752.31424518)
#define WGS84_E2    (0.00669437999014132)
#define OMEGA_E     (7.2921151467e-5)
#define mu (3.986005e14)
#define NUM_CHANS 10
#define PI (3.1415926535898)

const double SPEED_OF_LIGHT = 299792458.0;

int main(){
    double lat,lon,alt;
    struct Location predicted_location;
    struct Location *sv_location;
    int valid_locations = 0;
    int sv_ids[5]={1,2,3,4,5};
    int size = sizeof(sv_ids)/sizeof(int);

    sv_location = malloc(sizeof(struct Location)*size);
/*
space_vehicles[1].mean_motion_at_ephemeris = 0.0;
space_vehicles[1].sqrt_A = 0.515400631714e+04;
space_vehicles[1].Cus = 0.303611159325e-05;
space_vehicles[1].Cuc = 0.271014869213e-05;
space_vehicles[1].Cis = -0.275671482086e-06;
space_vehicles[1].Cic = -0.247731804848e-06;
space_vehicles[1].Crc = 0.337906250000e+03;
space_vehicles[1].Crs = 0.559687500000e+02;
space_vehicles[1].delta_n = 0.399659504565e-08;
space_vehicles[1].e =0.127465656260e-01;
space_vehicles[1].w = 0.0;
space_vehicles[1].omega_0 = -0.248667064980e+01;
space_vehicles[1].idot = 0.175007289756e-09;
space_vehicles[1].omega_dot = -0.796426031484e-08;
space_vehicles[1].inclination_at_ephemeris = 0.990713345450e+00;
space_vehicles[1].group_delay = 0.512227416039e-08;
space_vehicles[1].reference_time = 0.169008000000e+06;
space_vehicles[1].correction_f2 = 0.000000000000e+00;
space_vehicles[1].correction_f1 = 0.159161572810e-11;
space_vehicles[1].correction_f0 = 0.170257873833e-03;
space_vehicles[1].pos_t = 0.0;
space_vehicles[1].time_correction = 0.0;
space_vehicles[1].time_of_ephemeris = 0.172800000000e+06;
space_vehicles[1].lock_code_nco = 0.0;
space_vehicles[1].time_raw = 0.0;
space_vehicles[1].nav_subframe_week = 0.230200000000e+04;
space_vehicles[1].lock_ms_frame = 0.0;
space_vehicles[1].subframe ;
space_vehicles[1].subframe_number = 0.0;
space_vehicles[1].synced = 0.0;


// /*
space_vehicles[1].mean_motion_at_ephemeris = 0.0;
space_vehicles[1].sqrt_A = 0.0;
space_vehicles[1].Cus = 0.0;
space_vehicles[1].Cuc = 0.0;
space_vehicles[1].Cis = 0.0;
space_vehicles[1].Cic = 0.0;
space_vehicles[1].Crc = 0.0;
space_vehicles[1].Crs = 0.559687500000D+02;
space_vehicles[1].delta_n = 0.0;
space_vehicles[1].e = 0.0;
space_vehicles[1].w = 0.0;
space_vehicles[1].omega_0 = 0.0;
space_vehicles[1].idot = 0.0;
space_vehicles[1].omega_dot = 0.0;
space_vehicles[1].inclination_at_ephemeris = 0.0;
space_vehicles[1].group_delay = 0.0;
space_vehicles[1].reference_time = 0.0;
space_vehicles[1].correction_f2 = 0.0;
space_vehicles[1].correction_f1 = 0.0;
space_vehicles[1].correction_f0 = 0.0;
space_vehicles[1].pos_t = 0.0;
space_vehicles[1].time_correction = 0.0;
space_vehicles[1].time_of_ephemeris = 0.0;
space_vehicles[1].lock_code_nco = 0.0;
space_vehicles[1].time_raw = 0.0;
space_vehicles[1].nav_subframe_week = 0.0;
space_vehicles[1].lock_ms_frame = 0.0;
space_vehicles[1].subframe ;
space_vehicles[1].subframe_number = 0.0;
*/
    for (int i=0;i<size;i++){
        struct Space_vehicle *sv;
        sv = &space_vehicles[sv_ids[i]];
        unsigned phase_gold_code = (sv->lock_code_nco >>1) & 0x7FFFFFFF;
        sv->time_raw = sv->nav_subframe_week * 6000 + sv->lock_ms_frame + ((double)phase_gold_code)/(1023*(1<<21));
        sv->time_raw /= 1000.0;
        // nav_subframe(sv);
        sv_calc_corrected_time(i);
        sv_calc_location(i, sv_location + valid_locations);
        if(sv_location[valid_locations].x < 40000000 && sv_location[valid_locations].x > -40000000) 
            if(sv_location[valid_locations].y < 40000000 && sv_location[valid_locations].y > -40000000)
                if(sv_location[valid_locations].z < 40000000 && sv_location[valid_locations].z > -40000000)
                    valid_locations++;

    }
    if(valid_locations < 4){
        free(sv_location);
    }

    A_solve(valid_locations, sv_location, &predicted_location);
    LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z, lat, lon, alt);
    printf("Lat/Lon/Alt : %20.6f, %20.6f, %20.1f\n", lat*180/PI, lon*180/PI, alt);

    free(sv_location);
    return 0;
}



// Function definitions used above

void sv_calc_corrected_time(int i) {
    double delta_t, delta_tr, ek, time_correction;
    struct Space_vehicle *sv;    
    sv = space_vehicles+i;

    delta_t = sv->time_raw - sv->reference_time;
    if(delta_t >  302400)  delta_t -= 604800;
    if(delta_t < -302400)  delta_t += 604800;

    ek = orbit_ecc_anom(sv, sv->time_raw);

    delta_tr = -4.442807633e-10 * sv->e * sv->sqrt_A * sin(ek);

    time_correction = sv->correction_f0 
                    + (sv->correction_f1 * delta_t) 
                    + (sv->correction_f2 * delta_t * delta_t) 
                    + delta_tr
                    - sv->group_delay;
    sv->pos_t = sv->time_raw - time_correction; 

}

double  orbit_ecc_anom(struct Space_vehicle *sv, double t) {
    int    iterations  = 200;
    double delta       = pow(10,-10);
    double estimate, correction, semi_major_axis,     computed_mean_motion;
    double time_from_ephemeris,  correct_mean_motion, mean_anomaly;

    semi_major_axis      = sv->sqrt_A * sv->sqrt_A;
    computed_mean_motion = sqrt(mu / pow(semi_major_axis,3.0));

    time_from_ephemeris  = t - sv->time_of_ephemeris;
    if(time_from_ephemeris >  302400) time_from_ephemeris  -= 604800;
    if(time_from_ephemeris < -302400) time_from_ephemeris  += 604800;
    correct_mean_motion  = computed_mean_motion + sv->delta_n;

    /* Add on how much we have moved through the orbit since ephemeris */
    mean_anomaly         = sv->mean_motion_at_ephemeris + correct_mean_motion * time_from_ephemeris;

    /* First estimate */
    estimate   = (sv->e<0.8) ? mean_anomaly :  PI;
    correction = estimate - (mean_anomaly + sv->e*sin(mean_anomaly));

    /* Solve iteratively */
    while ((fabs(correction)>delta) && iterations > 0) {
        double last = estimate;
        estimate  = mean_anomaly  + sv->e * sin(estimate);
        correction = estimate - last;
        iterations--;
    }
    return estimate;
}

int sv_calc_location(int id, struct Location *l)
{  
    l->time = space_vehicles[id].pos_t;
    orbit_calc_position(space_vehicles+id, l);
    return 1;
}

int orbit_calc_position(struct Space_vehicle *sv, struct Location *l)
{
  double time_from_ephemeris,   semi_major_axis;
  double ek, true_anomaly,      corrected_argument_of_latitude;
  double argument_of_latitude,  argument_of_latitude_correction;
  double radius_correction,     corrected_radius;
  double correction_of_inclination;
  double pos_in_orbial_plane_x, pos_in_orbial_plane_y;
  double corrected_inclination, corrected_angle_of_ascending_node;
  
  /***********************
  * Calculate orbit
  ***********************/
  time_from_ephemeris  = sv->pos_t - sv->time_of_ephemeris;
  if(time_from_ephemeris >  302400) time_from_ephemeris  -= 604800;
  if(time_from_ephemeris < -302400) time_from_ephemeris  += 604800;

  semi_major_axis      = sv->sqrt_A * sv->sqrt_A;
  ek = orbit_ecc_anom(sv, sv->pos_t);

  /* Now calculate the first approximation of the latitude */
  true_anomaly = atan2( sqrt(1-sv->e * sv->e) * sin(ek), cos(ek) - sv->e);
  argument_of_latitude = true_anomaly + sv->w;

  /*****************************************
  * Second Harmonic Perbturbations 
  *****************************************/
  argument_of_latitude_correction = sv->Cus * sin(2*argument_of_latitude) 
                                  + sv->Cuc * cos(2*argument_of_latitude);

  radius_correction               = sv->Crc * cos(2*argument_of_latitude) 
                                  + sv->Crs * sin(2*argument_of_latitude);
  
  correction_of_inclination       = sv->Cic * cos(2*argument_of_latitude) 
                                  + sv->Cis * sin(2*argument_of_latitude);

  corrected_argument_of_latitude  = argument_of_latitude + argument_of_latitude_correction;
  corrected_radius                = semi_major_axis * (1- sv->e * cos(ek)) + radius_correction;
  corrected_inclination           = sv->inclination_at_ephemeris + correction_of_inclination 
                                  + sv->idot*time_from_ephemeris;

  pos_in_orbial_plane_x = corrected_radius * cos(corrected_argument_of_latitude);
  pos_in_orbial_plane_y = corrected_radius * sin(corrected_argument_of_latitude);
  

  corrected_angle_of_ascending_node = sv->omega_0
                                    + (sv->omega_dot - OMEGA_E)*time_from_ephemeris 
                                    - OMEGA_E * sv->time_of_ephemeris;

  /******************************************************
  * Project into Earth Centered, Earth Fixed coordinates
  ******************************************************/
  l->x = pos_in_orbial_plane_x * cos(corrected_angle_of_ascending_node)
       - pos_in_orbial_plane_y * cos(corrected_inclination) * sin(corrected_angle_of_ascending_node);
  l->y = pos_in_orbial_plane_x * sin(corrected_angle_of_ascending_node) 
       + pos_in_orbial_plane_y * cos(corrected_inclination) * cos(corrected_angle_of_ascending_node);
  l->z = pos_in_orbial_plane_y * sin(corrected_inclination);
  l->time = sv->pos_t;

  return 1;
}

int A_solve(int chans, struct Location *sv_l, struct Location *p) {
    int i, j, r, c;

    double t_tx[NUM_CHANS]; // Clock replicas in seconds since start of week

    double x_sv[NUM_CHANS],
           y_sv[NUM_CHANS],
           z_sv[NUM_CHANS];
  
    double t_pc;  // Uncorrected system time when clock replica snapshots taken
    double t_rx;    // Corrected GPS time

    double dPR[NUM_CHANS]; // Pseudo range error

    double jac[NUM_CHANS][4], ma[4][4], mb[4][4], mc[4][NUM_CHANS], md[4];

    double weight[NUM_CHANS];

    p->x = p->y = p->z = p->time = t_pc = 0.0;

    for (i=0; i<chans && i < NUM_CHANS; i++) {
        weight[i] = 1.0;
        x_sv[i]   = sv_l[i].x;
        y_sv[i]   = sv_l[i].y;
        z_sv[i]   = sv_l[i].z;
        t_tx[i]   = sv_l[i].time;
        t_pc     += sv_l[i].time;
    }

    // Approximate starting value for receiver clock
    t_pc = t_pc/chans + 75e-3;

    // Iterate to user xyzt solution using Taylor Series expansion:
  double err_mag;
    for(j=0; j<MAX_ITER; j++) {
        t_rx = t_pc - p->time;
        for (i=0; i<chans; i++) {
            // Convert SV position to ECI coords (20.3.3.4.3.3.2)
            double theta = (t_tx[i] - t_rx) * OMEGA_E;

            double x_sv_eci = x_sv[i]*cos(theta) - y_sv[i]*sin(theta);
            double y_sv_eci = x_sv[i]*sin(theta) + y_sv[i]*cos(theta);
            double z_sv_eci = z_sv[i];

            // Geometric range (20.3.3.4.3.4)
            double gr = sqrt(pow(p->x - x_sv_eci, 2) +
                             pow(p->y - y_sv_eci, 2) +
                             pow(p->z - z_sv_eci, 2));

            dPR[i] = SPEED_OF_LIGHT*(t_rx - t_tx[i]) - gr;

            jac[i][0] = (p->x - x_sv_eci) / gr;
            jac[i][1] = (p->y - y_sv_eci) / gr;
            jac[i][2] = (p->z - z_sv_eci) / gr;
            jac[i][3] = SPEED_OF_LIGHT;
        }

        // ma = transpose(H) * W * H
        for (r=0; r<4; r++)
            for (c=0; c<4; c++) {
            ma[r][c] = 0;
            for (i=0; i<chans; i++) ma[r][c] += jac[i][r]*weight[i]*jac[i][c];
        }

        double determinant =
            ma[0][3]*ma[1][2]*ma[2][1]*ma[3][0] - ma[0][2]*ma[1][3]*ma[2][1]*ma[3][0] - ma[0][3]*ma[1][1]*ma[2][2]*ma[3][0] + ma[0][1]*ma[1][3]*ma[2][2]*ma[3][0]+
            ma[0][2]*ma[1][1]*ma[2][3]*ma[3][0] - ma[0][1]*ma[1][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[1][2]*ma[2][0]*ma[3][1] + ma[0][2]*ma[1][3]*ma[2][0]*ma[3][1]+
            ma[0][3]*ma[1][0]*ma[2][2]*ma[3][1] - ma[0][0]*ma[1][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[1][0]*ma[2][3]*ma[3][1] + ma[0][0]*ma[1][2]*ma[2][3]*ma[3][1]+
            ma[0][3]*ma[1][1]*ma[2][0]*ma[3][2] - ma[0][1]*ma[1][3]*ma[2][0]*ma[3][2] - ma[0][3]*ma[1][0]*ma[2][1]*ma[3][2] + ma[0][0]*ma[1][3]*ma[2][1]*ma[3][2]+
            ma[0][1]*ma[1][0]*ma[2][3]*ma[3][2] - ma[0][0]*ma[1][1]*ma[2][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[2][0]*ma[3][3] + ma[0][1]*ma[1][2]*ma[2][0]*ma[3][3]+
            ma[0][2]*ma[1][0]*ma[2][1]*ma[3][3] - ma[0][0]*ma[1][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[1][0]*ma[2][2]*ma[3][3] + ma[0][0]*ma[1][1]*ma[2][2]*ma[3][3];

        // mb = inverse(ma) = inverse(transpose(H)*W*H)
        mb[0][0] = (ma[1][2]*ma[2][3]*ma[3][1] - ma[1][3]*ma[2][2]*ma[3][1] + ma[1][3]*ma[2][1]*ma[3][2] - ma[1][1]*ma[2][3]*ma[3][2] - ma[1][2]*ma[2][1]*ma[3][3] + ma[1][1]*ma[2][2]*ma[3][3]) / determinant;
        mb[0][1] = (ma[0][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[2][3]*ma[3][1] - ma[0][3]*ma[2][1]*ma[3][2] + ma[0][1]*ma[2][3]*ma[3][2] + ma[0][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[2][2]*ma[3][3]) / determinant;
        mb[0][2] = (ma[0][2]*ma[1][3]*ma[3][1] - ma[0][3]*ma[1][2]*ma[3][1] + ma[0][3]*ma[1][1]*ma[3][2] - ma[0][1]*ma[1][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[3][3] + ma[0][1]*ma[1][2]*ma[3][3]) / determinant;
        mb[0][3] = (ma[0][3]*ma[1][2]*ma[2][1] - ma[0][2]*ma[1][3]*ma[2][1] - ma[0][3]*ma[1][1]*ma[2][2] + ma[0][1]*ma[1][3]*ma[2][2] + ma[0][2]*ma[1][1]*ma[2][3] - ma[0][1]*ma[1][2]*ma[2][3]) / determinant;
        mb[1][0] = (ma[1][3]*ma[2][2]*ma[3][0] - ma[1][2]*ma[2][3]*ma[3][0] - ma[1][3]*ma[2][0]*ma[3][2] + ma[1][0]*ma[2][3]*ma[3][2] + ma[1][2]*ma[2][0]*ma[3][3] - ma[1][0]*ma[2][2]*ma[3][3]) / determinant;
        mb[1][1] = (ma[0][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[2][2]*ma[3][0] + ma[0][3]*ma[2][0]*ma[3][2] - ma[0][0]*ma[2][3]*ma[3][2] - ma[0][2]*ma[2][0]*ma[3][3] + ma[0][0]*ma[2][2]*ma[3][3]) / determinant;
        mb[1][2] = (ma[0][3]*ma[1][2]*ma[3][0] - ma[0][2]*ma[1][3]*ma[3][0] - ma[0][3]*ma[1][0]*ma[3][2] + ma[0][0]*ma[1][3]*ma[3][2] + ma[0][2]*ma[1][0]*ma[3][3] - ma[0][0]*ma[1][2]*ma[3][3]) / determinant;
        mb[1][3] = (ma[0][2]*ma[1][3]*ma[2][0] - ma[0][3]*ma[1][2]*ma[2][0] + ma[0][3]*ma[1][0]*ma[2][2] - ma[0][0]*ma[1][3]*ma[2][2] - ma[0][2]*ma[1][0]*ma[2][3] + ma[0][0]*ma[1][2]*ma[2][3]) / determinant;
        mb[2][0] = (ma[1][1]*ma[2][3]*ma[3][0] - ma[1][3]*ma[2][1]*ma[3][0] + ma[1][3]*ma[2][0]*ma[3][1] - ma[1][0]*ma[2][3]*ma[3][1] - ma[1][1]*ma[2][0]*ma[3][3] + ma[1][0]*ma[2][1]*ma[3][3]) / determinant;
        mb[2][1] = (ma[0][3]*ma[2][1]*ma[3][0] - ma[0][1]*ma[2][3]*ma[3][0] - ma[0][3]*ma[2][0]*ma[3][1] + ma[0][0]*ma[2][3]*ma[3][1] + ma[0][1]*ma[2][0]*ma[3][3] - ma[0][0]*ma[2][1]*ma[3][3]) / determinant;
        mb[2][2] = (ma[0][1]*ma[1][3]*ma[3][0] - ma[0][3]*ma[1][1]*ma[3][0] + ma[0][3]*ma[1][0]*ma[3][1] - ma[0][0]*ma[1][3]*ma[3][1] - ma[0][1]*ma[1][0]*ma[3][3] + ma[0][0]*ma[1][1]*ma[3][3]) / determinant;
        mb[2][3] = (ma[0][3]*ma[1][1]*ma[2][0] - ma[0][1]*ma[1][3]*ma[2][0] - ma[0][3]*ma[1][0]*ma[2][1] + ma[0][0]*ma[1][3]*ma[2][1] + ma[0][1]*ma[1][0]*ma[2][3] - ma[0][0]*ma[1][1]*ma[2][3]) / determinant;
        mb[3][0] = (ma[1][2]*ma[2][1]*ma[3][0] - ma[1][1]*ma[2][2]*ma[3][0] - ma[1][2]*ma[2][0]*ma[3][1] + ma[1][0]*ma[2][2]*ma[3][1] + ma[1][1]*ma[2][0]*ma[3][2] - ma[1][0]*ma[2][1]*ma[3][2]) / determinant;
        mb[3][1] = (ma[0][1]*ma[2][2]*ma[3][0] - ma[0][2]*ma[2][1]*ma[3][0] + ma[0][2]*ma[2][0]*ma[3][1] - ma[0][0]*ma[2][2]*ma[3][1] - ma[0][1]*ma[2][0]*ma[3][2] + ma[0][0]*ma[2][1]*ma[3][2]) / determinant;
        mb[3][2] = (ma[0][2]*ma[1][1]*ma[3][0] - ma[0][1]*ma[1][2]*ma[3][0] - ma[0][2]*ma[1][0]*ma[3][1] + ma[0][0]*ma[1][2]*ma[3][1] + ma[0][1]*ma[1][0]*ma[3][2] - ma[0][0]*ma[1][1]*ma[3][2]) / determinant;
        mb[3][3] = (ma[0][1]*ma[1][2]*ma[2][0] - ma[0][2]*ma[1][1]*ma[2][0] + ma[0][2]*ma[1][0]*ma[2][1] - ma[0][0]*ma[1][2]*ma[2][1] - ma[0][1]*ma[1][0]*ma[2][2] + ma[0][0]*ma[1][1]*ma[2][2]) / determinant;

        // mc = inverse(transpose(H)*W*H) * transpose(H)
        for (r=0; r<4; r++)
            for (c=0; c<chans; c++) {
            mc[r][c] = 0;
            for (i=0; i<4; i++) mc[r][c] += mb[r][i]*jac[c][i];
        }

        // md = inverse(transpose(H)*W*H) * transpose(H) * W * dPR
        for (r=0; r<4; r++) {
            md[r] = 0;
            for (i=0; i<chans; i++) md[r] += mc[r][i]*weight[i]*dPR[i];
        }

        double dx = md[0];
        double dy = md[1];
        double dz = md[2];
        double dt = md[3];

        err_mag = sqrt(dx*dx + dy*dy + dz*dz);


        if (err_mag<1.0) break;

        p->x    += dx;
        p->y    += dy;
        p->z    += dz;
        p->time += dt;
    }
    return j;
}

 void LatLonAlt(
    double x_n, double y_n, double z_n,
    double lat, double lon, double alt) {
    int iterations = 100;
    const double a  = WGS84_A;
    const double e2 = WGS84_E2;

    const double p = sqrt(x_n*x_n + y_n*y_n);

    lon = 2.0 * atan2(y_n, x_n + p);
    lat = atan(z_n / (p * (1.0 - e2)));
    alt = 0.0;

    while(iterations > 0) {
        double tmp = alt;
        double N = a / sqrt(1.0 - e2*pow(sin(lat),2));
        alt = p/cos(lat) - N;
        lat = atan(z_n / (p * (1.0 - e2*N/(N + alt))));
        if(fabs(alt-tmp)<1e-3) 
            break;
        iterations--;
    }
}
