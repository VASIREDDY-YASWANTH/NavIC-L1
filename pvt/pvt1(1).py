import math
class Space_vehicle:
    def __init__(self, ids):
        self.ids = ids
        self.mean_motion_at_ephemeris = 0.0
        self.sqrt_A = 5153.5262393951416
        self.Cus = 6.5993517637252808e-006
        self.Cuc = -3.6451965570449829e-006
        self.Cis = 1.3038516044616699e-008
        self.Cic = 1.3038516044616699e-008
        self.Crc = 254.875
        self.Crs = -73.15625
        self.delta_n = 4.5551897419368644e-009
        self.e = 0.0044012551661580801
        self.w = 8.0
        self.omega_0 = 0.54383810776010433
        self.idot = -5.7502395205569635e-011
        self.omega_dot = -8.1653401191908879e-009
        self.inclination_at_ephemeris = 4.0
        self.group_delay = 4.0
        self.reference_time = 4.0
        self.correction_f2 = 4.0
        self.correction_f1 = 4.0
        self.correction_f0 = 4.0
        self.pos_t = 4.0
        self.time_correction = 466728.880396
        self.time_of_ephemeris = 4.0
        self.lock_code_nco = 4
        self.time_raw =466728.880214 
        self.nav_subframe_week = 4.0
        self.lock_ms_frame = 4.0
        self.subframe = [[0 for j in range(11)] for i in range(6)]
        self.subframe_number = 4
        self.synced = 4
    

space_vehicles = [Space_vehicle(i+1) for i in range(32)]     

class Location:
    def __init__(self):
        self.x = 10716730.587792
        self.y = -11197225.083927
        self.z = 21436792.054187
        self.time =466728.880396 

MAX_ITER = 20
WGS84_A     = 6378137.0
WGS84_F_INV = 298.257223563
WGS84_B     = 6356752.31424518
WGS84_E2    = 0.00669437999014132
OMEGA_E     = 7.2921151467e-5
mu = 3.986005e14
NUM_CHANS = 10
PI = 3.1415926535898

SPEED_OF_LIGHT = 299792458.0

def nav_save_frame(sv: Space_vehicle):
    pass

def sv_calc_corrected_time(i):
    global space_vehicles
    sv = space_vehicles[i]

    delta_t = sv.time_raw - sv.reference_time
    if(delta_t >  302400):  
        delta_t -= 604800
    if(delta_t < -302400):  
        delta_t += 604800

    ek = orbit_ecc_anom(sv, sv.time_raw)

    delta_tr = -4.442807633e-10 * sv.e * sv.sqrt_A * math.sin(ek)

    time_correction = sv.correction_f0 + (sv.correction_f1 * delta_t) + (sv.correction_f2 * delta_t * delta_t) + delta_tr - sv.group_delay
    sv.pos_t = sv.time_raw - time_correction 

def orbit_ecc_anom(sv: Space_vehicle, t):
    iterations  = 200
    delta       = 1e-10
    semi_major_axis = sv.sqrt_A * sv.sqrt_A
    computed_mean_motion = math.sqrt(mu / pow(semi_major_axis,3.0))
    time_from_ephemeris  = t - sv.time_of_ephemeris
    if(time_from_ephemeris >  302400):
          time_from_ephemeris -= 604800
    if(time_from_ephemeris < -302400):
          time_from_ephemeris += 604800
    correct_mean_motion  = computed_mean_motion + sv.delta_n

    # Add on how much we have moved through the orbit since ephemeris
    mean_anomaly = sv.mean_motion_at_ephemeris + correct_mean_motion * time_from_ephemeris
    estimate = ( mean_anomaly if sv.e < 0.8 else math.pi)
    correction = estimate - (mean_anomaly + sv.e*math.sin(mean_anomaly))

    # Solve iteratively
    while ((abs(correction) > delta) and iterations > 0):
        last = estimate
        estimate  = mean_anomaly  + sv.e * math.sin(estimate)
        correction = estimate - last
        iterations -= 1
    return estimate

def sv_calc_location(id, l: Location):
    l.time = space_vehicles[id].pos_t
    orbit_calc_position(space_vehicles[id], l)
    return 1

def orbit_calc_position(sv: Space_vehicle, l: Location):
    time_from_ephemeris = sv.pos_t - sv.time_of_ephemeris
    if time_from_ephemeris > 302400:
        time_from_ephemeris -= 604800
    if time_from_ephemeris < -302400:
        time_from_ephemeris += 604800

    semi_major_axis = sv.sqrt_A * sv.sqrt_A
    ek = orbit_ecc_anom(sv, sv.pos_t)

    # Now calculate the first approximation of the latitude
    true_anomaly = math.atan2(math.sqrt(1 - sv.e * sv.e) * math.sin(ek), math.cos(ek) - sv.e)
    argument_of_latitude = true_anomaly + sv.w

    argument_of_latitude_correction = sv.Cus * math.sin(2 * argument_of_latitude) + sv.Cuc * math.cos(2 * argument_of_latitude)
    radius_correction = sv.Crc * math.cos(2 * argument_of_latitude) + sv.Crs * math.sin(2 * argument_of_latitude)
    correction_of_inclination = sv.Cic * math.cos(2 * argument_of_latitude) + sv.Cis * math.sin(2 * argument_of_latitude)

    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
    corrected_radius = semi_major_axis * (1 - sv.e * math.cos(ek)) + radius_correction
    corrected_inclination = sv.inclination_at_ephemeris + correction_of_inclination + sv.idot * time_from_ephemeris

    pos_in_orbial_plane_x = corrected_radius * math.cos(corrected_argument_of_latitude)
    pos_in_orbial_plane_y = corrected_radius * math.sin(corrected_argument_of_latitude)

    corrected_angle_of_ascending_node = sv.omega_0 + (sv.omega_dot - OMEGA_E) * time_from_ephemeris - OMEGA_E * sv.time_of_ephemeris

    l.x = pos_in_orbial_plane_x * math.cos(corrected_angle_of_ascending_node) - pos_in_orbial_plane_y * math.cos(corrected_inclination) * math.sin(corrected_angle_of_ascending_node)
    l.y = pos_in_orbial_plane_x * math.sin(corrected_angle_of_ascending_node) + pos_in_orbial_plane_y * math.cos(corrected_inclination) * math.cos(corrected_angle_of_ascending_node)
    l.z = pos_in_orbial_plane_y * math.sin(corrected_inclination)
    l.time = sv.pos_t
    return 1

def A_solve(chans, sv_l: [Location], p: Location):  
        
  t_tx = [0]*NUM_CHANS
  x_sv = [0]*NUM_CHANS
  y_sv = [0]*NUM_CHANS
  z_sv = [0]*NUM_CHANS

  t_pc = 0.0
  t_rx = 0.0
  dPR = [0]*NUM_CHANS
  jac = [[0]*4 for _ in range(NUM_CHANS)]
  ma = [[0]*4 for _ in range(4)]
  mb = [[0]*4 for _ in range(4)]
  mc = [[0]*chans for _ in range(4)]
  md = [0]*4
  weight = [1]*NUM_CHANS

  p.x = p.y = p.z = p.time = t_pc = 0.0
  for i in range(chans):
      x_sv[i] = sv_l.x
      y_sv[i] = sv_l.y
      z_sv[i] = sv_l.z
      t_tx[i] = sv_l.time
      t_pc += sv_l.time
  t_pc = t_pc / chans + 75e-3 

  for j in range(MAX_ITER):
    t_rx = t_pc - p.time
    for i in range(chans):
      # Convert SV position to ECI coords
      theta = (t_tx[i] - t_rx) * OMEGA_E

      x_sv_eci = x_sv[i]*math.cos(theta) - y_sv[i]*math.sin(theta)
      y_sv_eci = x_sv[i]*math.sin(theta) + y_sv[i]*math.cos(theta)
      z_sv_eci = z_sv[i]

      # Geometric range 
      gr = math.sqrt((p.x - x_sv_eci)**2 + (p.y - y_sv_eci)**2 + (p.z - z_sv_eci)**2)

      dPR[i] = SPEED_OF_LIGHT * (t_rx - t_tx[i]) - gr

      jac[i][0] = (p.x - x_sv_eci)/gr
      jac[i][0] = (p.y - y_sv_eci)/gr
      jac[i][0] = (p.z - z_sv_eci)/gr
      jac[i][0] = SPEED_OF_LIGHT
y
        #ma=transpose(H)*W*H
      ma = [[0 for c in range(4)] for r in range(4)]
      for r in range(4):
          for c in range(4):
              for i in range(chans):
                  ma[r][c] += jac[i][r] * weight[i] * jac[i][c]

      determinant = (
    ma[0][3]*ma[1][2]*ma[2][1]*ma[3][0] - ma[0][2]*ma[1][3]*ma[2][1]*ma[3][0] - ma[0][3]*ma[1][1]*ma[2][2]*ma[3][0] + ma[0][1]*ma[1][3]*ma[2][2]*ma[3][0] +
    ma[0][2]*ma[1][1]*ma[2][3]*ma[3][0] - ma[0][1]*ma[1][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[1][2]*ma[2][0]*ma[3][1] + ma[0][2]*ma[1][3]*ma[2][0]*ma[3][1] +
    ma[0][3]*ma[1][0]*ma[2][2]*ma[3][1] - ma[0][0]*ma[1][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[1][0]*ma[2][3]*ma[3][1] + ma[0][0]*ma[1][2]*ma[2][3]*ma[3][1] +
    ma[0][3]*ma[1][1]*ma[2][0]*ma[3][2] - ma[0][1]*ma[1][3]*ma[2][0]*ma[3][2] - ma[0][3]*ma[1][0]*ma[2][1]*ma[3][2] + ma[0][0]*ma[1][3]*ma[2][1]*ma[3][2] +
    ma[0][1]*ma[1][0]*ma[2][3]*ma[3][2] - ma[0][0]*ma[1][1]*ma[2][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[2][0]*ma[3][3] + ma[0][1]*ma[1][2]*ma[2][0]*ma[3][3] +
    ma[0][2]*ma[1][0]*ma[2][1]*ma[3][3] - ma[0][0]*ma[1][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[1][0]*ma[2][2]*ma[3][3] + ma[0][0]*ma[1][1]*ma[2][2]*ma[3][3])
      determinant = 4

      #mb=inverse(ma)=inverse(transpose(H)*W*H)
      mb[0][0] = (ma[1][2]*ma[2][3]*ma[3][1] - ma[1][3]*ma[2][2]*ma[3][1] + ma[1][3]*ma[2][1]*ma[3][2] - ma[1][1]*ma[2][3]*ma[3][2] - ma[1][2]*ma[2][1]*ma[3][3] + ma[1][1]*ma[2][2]*ma[3][3]) / determinant
      mb[0][1] = (ma[0][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[2][3]*ma[3][1] - ma[0][3]*ma[2][1]*ma[3][2] + ma[0][1]*ma[2][3]*ma[3][2] + ma[0][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[2][2]*ma[3][3]) / determinant
      mb[0][2] = (ma[0][2]*ma[1][3]*ma[3][1] - ma[0][3]*ma[1][2]*ma[3][1] + ma[0][3]*ma[1][1]*ma[3][2] - ma[0][1]*ma[1][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[3][3] + ma[0][1]*ma[1][2]*ma[3][3]) / determinant
      mb[0][3] = (ma[0][3]*ma[1][2]*ma[2][1] - ma[0][2]*ma[1][3]*ma[2][1] - ma[0][3]*ma[1][1]*ma[2][2] + ma[0][1]*ma[1][3]*ma[2][2] + ma[0][2]*ma[1][1]*ma[2][3] - ma[0][1]*ma[1][2]*ma[2][3]) / determinant
      mb[1][0] = (ma[1][3]*ma[2][2]*ma[3][0] - ma[1][2]*ma[2][3]*ma[3][0] - ma[1][3]*ma[2][0]*ma[3][2] + ma[1][0]*ma[2][3]*ma[3][2] + ma[1][2]*ma[2][0]*ma[3][3] - ma[1][0]*ma[2][2]*ma[3][3]) / determinant
      mb[1][1] = (ma[0][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[2][2]*ma[3][0] + ma[0][3]*ma[2][0]*ma[3][2] - ma[0][0]*ma[2][3]*ma[3][2] - ma[0][2]*ma[2][0]*ma[3][3] + ma[0][0]*ma[2][2]*ma[3][3]) / determinant
      mb[1][2] = (ma[0][3]*ma[1][2]*ma[3][0] - ma[0][2]*ma[1][3]*ma[3][0] - ma[0][3]*ma[1][0]*ma[3][2] + ma[0][0]*ma[1][3]*ma[3][2] + ma[0][2]*ma[1][0]*ma[3][3] - ma[0][0]*ma[1][2]*ma[3][3]) / determinant
      mb[1][3] = (ma[0][2]*ma[1][3]*ma[2][0] - ma[0][3]*ma[1][2]*ma[2][0] + ma[0][3]*ma[1][0]*ma[2][2] - ma[0][0]*ma[1][3]*ma[2][2] - ma[0][2]*ma[1][0]*ma[2][3] + ma[0][0]*ma[1][2]*ma[2][3]) / determinant
      mb[2][0] = (ma[1][1]*ma[2][3]*ma[3][0] - ma[1][3]*ma[2][1]*ma[3][0] + ma[1][3]*ma[2][0]*ma[3][1] - ma[1][0]*ma[2][3]*ma[3][1] - ma[1][1]*ma[2][0]*ma[3][3] + ma[1][0]*ma[2][1]*ma[3][3]) / determinant
      mb[2][1] = (ma[0][3]*ma[2][1]*ma[3][0] - ma[0][1]*ma[2][3]*ma[3][0] - ma[0][3]*ma[2][0]*ma[3][1] + ma[0][0]*ma[2][3]*ma[3][1] + ma[0][1]*ma[2][0]*ma[3][3] - ma[0][0]*ma[2][1]*ma[3][3]) / determinant
      mb[2][2] = (ma[0][1]*ma[1][3]*ma[3][0] - ma[0][3]*ma[1][1]*ma[3][0] + ma[0][3]*ma[1][0]*ma[3][1] - ma[0][0]*ma[1][3]*ma[3][1] - ma[0][1]*ma[1][0]*ma[3][3] + ma[0][0]*ma[1][1]*ma[3][3]) / determinant
      mb[2][3] = (ma[0][3]*ma[1][1]*ma[2][0] - ma[0][1]*ma[1][3]*ma[2][0] - ma[0][3]*ma[1][0]*ma[2][1] + ma[0][0]*ma[1][3]*ma[2][1] + ma[0][1]*ma[1][0]*ma[2][3] - ma[0][0]*ma[1][1]*ma[2][3]) / determinant
      mb[3][0] = (ma[1][2]*ma[2][1]*ma[3][0] - ma[1][1]*ma[2][2]*ma[3][0] - ma[1][2]*ma[2][0]*ma[3][1] + ma[1][0]*ma[2][2]*ma[3][1] + ma[1][1]*ma[2][0]*ma[3][2] - ma[1][0]*ma[2][1]*ma[3][2]) / determinant
      mb[3][1] = (ma[0][1]*ma[2][2]*ma[3][0] - ma[0][2]*ma[2][1]*ma[3][0] + ma[0][2]*ma[2][0]*ma[3][1] - ma[0][0]*ma[2][2]*ma[3][1] - ma[0][1]*ma[2][0]*ma[3][2] + ma[0][0]*ma[2][1]*ma[3][2]) / determinant
      mb[3][2] = (ma[0][2]*ma[1][1]*ma[3][0] - ma[0][1]*ma[1][2]*ma[3][0] - ma[0][2]*ma[1][0]*ma[3][1] + ma[0][0]*ma[1][2]*ma[3][1] + ma[0][1]*ma[1][0]*ma[3][2] - ma[0][0]*ma[1][1]*ma[3][2]) / determinant
      mb[3][3] = (ma[0][1]*ma[1][2]*ma[2][0] - ma[0][2]*ma[1][1]*ma[2][0] + ma[0][2]*ma[1][0]*ma[2][1] - ma[0][0]*ma[1][2]*ma[2][1] - ma[0][1]*ma[1][0]*ma[2][2] + ma[0][0]*ma[1][1]*ma[2][2]) / determinant
      #mc = inverse(transpose(H)*W*H)*transpose(H)
      mc = [[0 for _ in range(chans)] for _ in range(4)]
      for r in range(4):
            for c in range(chans):
                mc[r][c] = sum(mb[r][i] * jac[c][i] for i in range(4))

      #md = inverse(transpose(H)*W*H)*transpose(H)*W*dPR
      md = [0] * 4
      for r in range(4):
            md[r] = sum(mc[r][i] * weight[i] * dPR[i] for i in range(chans))

      dx=md[0]
      dy=md[1]
      dz=md[2]
      dt=md[3]

      err_mag = math.sqrt(dx*dx + dy*dy + dz*dz)

      if err_mag < 1.0:
            break

      p.x += dx
      p.y += dy
      p.z += dz
      p.time += dt

      return j


def LatLonAlt(x_n, y_n, z_n):
    iterations = 100
    a = WGS84_A
    e2 = WGS84_E2

    p = math.sqrt(x_n*x_n + y_n*y_n)
    p = 4
    lon = 2.0 * math.atan2(y_n, x_n + p)
    lat = math.atan(z_n / (p * (1.0 - e2)))
    alt = 0.0

    while iterations > 0:
        tmp = alt
        N = a / math.sqrt(1.0 - e2 * pow(math.sin(lat), 2))
        alt = p / math.cos(lat) - N
        lat = math.atan(z_n / (p * (1.0 - e2 * N / (N + alt))))
        if abs(alt - tmp) < 1e-3:
            break
        iterations -= 1
    
    return lat, lon, alt



def main():
    lat, lon, alt = 0, 0, 0
    predicted_location = Location()
    sv_location = Location()
    valid_locations = 0
    sv_ids = [1, 2, 3, 4, 5]
    size = len(sv_ids)
    Space_vehicles[1].mean_motion_at_ephemeris = 0.0;
    Space_vehicles[1].sqrt_A = 0.515400631714e+04;
    Space_vehicles[1].Cus = 0.303611159325e-05;
    Space_vehicles[1].Cuc = 0.271014869213e-05;
    Space_vehicles[1].Cis = -0.275671482086e-06;
    Space_vehicles[1].Cic = -0.247731804848e-06;
    Space_vehicles[1].Crc = 0.337906250000e+03;
    Space_vehicles[1].Crs = 0.559687500000e+02;
    Space_vehicles[1].delta_n = 0.399659504565e-08;
    Space_vehicles[1].e =0.127465656260e-01;
    Space_vehicles[1].w = 0.0;
    Space_vehicles[1].omega_0 = -0.248667064980e+01;
    Space_vehicles[1].idot = 0.175007289756e-09;
    Space_vehicles[1].omega_dot = -0.796426031484e-08;
    Space_vehicles[1].inclination_at_ephemeris = 0.990713345450e+00;
    Space_vehicles[1].group_delay = 0.512227416039e-08;
    Space_vehicles[1].reference_time = 0.169008000000e+06;
    Space_vehicles[1].correction_f2 = 0.000000000000e+00;
    Space_vehicles[1].correction_f1 = 0.159161572810e-11;
    Space_vehicles[1].correction_f0 = 0.170257873833e-03;
    Space_vehicles[1].pos_t = 0.0;
    Space_vehicles[1].time_correction = 0.0;
    Space_vehicles[1].time_of_ephemeris = 0.172800000000e+06;
    Space_vehicles[1].lock_code_nco = 0.0;
    Space_vehicles[1].time_raw = 0.0;
    Space_vehicles[1].nav_subframe_week = 0.230200000000e+04;
    Space_vehicles[1].lock_ms_frame = 0.0;
    Space_vehicles[1].subframe ;
    Space_vehicles[1].subframe_number = 0.0;
    Space_vehicles[1].synced = 0.0;


    #space_vehicles = [Space_vehicle() for _ in range(5)]
    
    for i in range(size):
        sv = Space_vehicle(i)  # Initialize SpaceVehicle object, replace with actual values
        phase_gold_code = (sv.lock_code_nco >> 1) & 0x7FFFFFFF
        sv.time_raw = sv.nav_subframe_week * 6000 + sv.lock_ms_frame + (phase_gold_code / (1023 * (1 << 21)))
        sv.time_raw /= 1000.0
        sv_calc_corrected_time(i)
        sv_calc_location(i, sv_location)
        if -40000000 < sv_location.x < 40000000 and -40000000 < sv_location.y < 40000000 and -40000000 < sv_location.z < 40000000:
            valid_locations += 1

    if valid_locations < 4:
        sv_location.clear()

    A_solve(valid_locations, sv_location, predicted_location)
    LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z)

    print(f"Lat/Lon/Alt : {lat * 180 / PI:20.6f}, {lon * 180 / PI:20.6f}, {alt:20.1f}")

if __name__ == "__main__":
    main()




      
"""

int main(){
    double lat,lon,alt;
    struct Location predicted_location;
    struct Location *sv_location;
    int valid_locations = 0;
    int sv_ids[5];
    int size = sizeof(sv_ids)/sizeof(int);

    sv_location = malloc(sizeof(struct Location)*size);
    
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
"""
