import math
#class defining the space vehicle parameters
class Space_vehicle:
    def __init__(self, ids,mmae,s_A,cus,cuc,cis,cic,crc,crs,d_n,e,w,omega_0,idot,omega_dot,iae,g_d,r_t,f2,f1,f0,pos_t,toe,lock_code_nco,time_raw,nav_subframe_week,lock_ms_frame,synced):
        self.ids = ids
        self.mean_motion_at_ephemeris = mmae
        self.sqrt_A = s_A
        self.Cus = cus
        self.Cuc = cuc
        self.Cis = cis
        self.Cic = cic
        self.Crc = crc
        self.Crs = crs
        self.delta_n = d_n
        self.e = e
        self.w = w
        self.omega_0 = omega_0
        self.idot = idot
        self.omega_dot = omega_dot
        self.inclination_at_ephemeris = iae
        self.group_delay = g_d
        self.reference_time = r_t
        self.correction_f2 = f2
        self.correction_f1 = f1
        self.correction_f0 = f0
        self.pos_t = pos_t
        self.time_of_ephemeris = toe
        self.lock_code_nco = lock_code_nco
        self.time_raw = time_raw
        self.nav_subframe_week = nav_subframe_week
        self.lock_ms_frame = lock_ms_frame
        self.synced = synced
#class defining the location parameters of each space vehicle
class Location:
    def __init__(self, X , Y , Z, T):
        self.x = X
        self.y = Y
        self.z = Z
        self.time =T
#defining CONSTANTS
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

#initializing space vehicle parameters
S_V=[]
#space vehicle 1 parameters
S_V.append(Space_vehicle(1, 0.627712271718553, 5153.52623939514, 0.00000659935176372528, -0.00000364519655704498, 0.00000001303851604462, 0.00000001303851604462, 254.875, -73.15625, 0.00000000455518974194, 0.00440125516615808, 0.948704618521089, 0.543838107760104, -0.00000000005750239521, -0.00000000816534011919, 0.96682775824546, -0.00000001909211277962, 468000, 0, -0.00000000000363797881, -0.000181447714567184, 0, 468000,    872457844  , 0  , 77786, 4, 0)) 
#space vehicle 2 parameters
S_V.append(Space_vehicle(4, 3.02458803990085, 5153.63453292847, 0.00000938959419727325, -0.00000134296715259552, 0.00000004842877388, 0.00000005215406417847, 196.9375, -25.53125, 0.00000000449018703431, 0.00251956866122782, -1.25409232099612, 1.62825900432321, 0.00000000013464846579,  -0.00000000786639909567, 0.960616162730431, -0.00000000838190317154, 468000, 0, 0.00000000000284217094, 0.0001630624756217, 0, 468000,   101717271  , 0,  77786 , 5,  0  )) 
#space vehicle 3 parameters
S_V.append(Space_vehicle(5, 0.796297202986197, 5155.34639358521, 0.0000064503401517868, -0.00000375322997570038, -0.00000006891787052155, -0.00000016018748283386, 256.9375, -70.65625, 0.00000000480912889095, 0.0130847115069628, 1.55579952931809, 0.503145066445881, -0.00000000005893102614, -0.00000000841285042899, 0.958335540739085, -0.00000000791624188423, 468000, 0, 0.00000000000341060513, 0.000295627396553755, 0, 468000,  3532408096   , 0,  77786     , 6,  0 )) 
#space vehicle 4 parameters
S_V.append(Space_vehicle(6, -2.60055087716526, 5153.80060768127, 0.00000423379242420197, 0.00000039488077163696, 0.00000013038516044617, -0.00000001490116119385, 305.5625, 3.8125, 0.00000000427982112887, 0.00755723693873733, -1.04896671410838, -0.494186142962186, 0.00000000016893560827, -0.00000000794640242813, 0.978481871297476, -0.00000001303851604462, 467984, 0, 0.00000000000375166564, 0.00000797677785158157, 0, 467984,  2183385124  , 0,   77786  , 3, 0))  

#default initilzation of location parameters of space vehicles to zero 
L=[]
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))

#************************************************************************
#********************************MAIN CODE*******************************
#************************************************************************
def main():
    lat, lon, alt = 0, 0, 0
    predicted_location = Location(0,0,0,0)
    sv_location = Location(0,0,0,0)
    valid_locations = 0
    sv_ids = [1, 2, 3, 4 ]
    size = len(sv_ids)
    
    for i in range(size):
        sv = S_V[i]  
        sv_location=L[i]
        phase_gold_code = (sv.lock_code_nco >> 1) & 0x7FFFFFFF
        #calculating transmit time in milliseconds
        sv.time_raw = sv.nav_subframe_week * 6000 + sv.lock_ms_frame + (float(phase_gold_code) / (1023 * (1 << 21)))
        #converting from milliseconds to seconds
        sv.time_raw /= 1000.0
        #correcting the time using calibration factors
        sv_calc_corrected_time(sv)
        #calculating the location of the space vehicle
        sv_calc_location(sv, sv_location)
        #verifying that the locations are valid or not
        if -40000000 < sv_location.x < 40000000 and -40000000 < sv_location.y < 40000000 and -40000000 < sv_location.z < 40000000:
            valid_locations += 1

        if i == 3:
         A_solve(valid_locations,L , predicted_location)
         lat,lon,alt = LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z)
         print(f"Lat/Lon/Alt : {lat * 180 / PI:20.6f}, {lon * 180 / PI:20.6f}, {alt:20.1f}")




#************************************************************************
def sv_calc_corrected_time(sv):
    #calculate the time for the adjustment
    delta_t = sv.time_raw - sv.reference_time
    if(delta_t >  302400):  
        delta_t -= 604800
    if(delta_t < -302400):  
        delta_t += 604800
    #Relativtstic term
    ek = orbit_ecc_anom(sv, sv.time_raw)
    delta_tr = -4.442807633e-10 * sv.e * sv.sqrt_A * math.sin(ek)
    time_correction = sv.correction_f0 + (sv.correction_f1 * delta_t) + (sv.correction_f2 * delta_t * delta_t) + delta_tr - sv.group_delay
    sv.pos_t = sv.time_raw - time_correction

#************************************************************************
def orbit_ecc_anom(sv, t):
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
    #Add on how much we have moved through the orbit since ephemeris
    mean_anomaly = sv.mean_motion_at_ephemeris + correct_mean_motion * time_from_ephemeris
    #First estimate
    estimate = ( mean_anomaly if sv.e < 0.8 else math.pi)
    correction = estimate - (mean_anomaly + sv.e*math.sin(mean_anomaly))
    #Solve iteratively
    while ((abs(correction) > delta) and iterations > 0):
        last = estimate
        estimate  = mean_anomaly  + sv.e * math.sin(estimate)
        correction = estimate - last
        iterations -= 1
    return estimate

#************************************************************************
def sv_calc_location(sv, l: Location):
    l.time = sv.pos_t
    orbit_calc_position(sv, l)
    return 1

#************************************************************************
#*******Calculate where the space vehicle will be at time "pos_t"********
def orbit_calc_position(sv, l: Location):
    #calculate orbit
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
    #second Harmonic Perbturbations
    argument_of_latitude_correction = sv.Cus * math.sin(2 * argument_of_latitude) + sv.Cuc * math.cos(2 * argument_of_latitude)
    radius_correction = sv.Crc * math.cos(2 * argument_of_latitude) + sv.Crs * math.sin(2 * argument_of_latitude)
    correction_of_inclination = sv.Cic * math.cos(2 * argument_of_latitude) + sv.Cis * math.sin(2 * argument_of_latitude)

    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
    corrected_radius = semi_major_axis * (1 - sv.e * math.cos(ek)) + radius_correction
    corrected_inclination = sv.inclination_at_ephemeris + correction_of_inclination + sv.idot * time_from_ephemeris

    pos_in_orbial_plane_x = corrected_radius * math.cos(corrected_argument_of_latitude)
    pos_in_orbial_plane_y = corrected_radius * math.sin(corrected_argument_of_latitude)

    corrected_angle_of_ascending_node = sv.omega_0 + (sv.omega_dot - OMEGA_E) * time_from_ephemeris - OMEGA_E * sv.time_of_ephemeris
    #project into earth centered ,earth fixed coordinates
    l.x = pos_in_orbial_plane_x * math.cos(corrected_angle_of_ascending_node) - pos_in_orbial_plane_y * math.cos(corrected_inclination) * math.sin(corrected_angle_of_ascending_node)
    l.y = pos_in_orbial_plane_x * math.sin(corrected_angle_of_ascending_node) + pos_in_orbial_plane_y * math.cos(corrected_inclination) * math.cos(corrected_angle_of_ascending_node)
    l.z = pos_in_orbial_plane_y * math.sin(corrected_inclination)
    l.time = sv.pos_t
    return 1

#************************************************************************
def A_solve(chans, sv_l , p ):  
  t_tx = [0]*NUM_CHANS #clock replicas in second since start of week
  x_sv = [0]*NUM_CHANS
  y_sv = [0]*NUM_CHANS
  z_sv = [0]*NUM_CHANS
  
  t_pc = 0.0 #uncorrected system time when clock replica snapshots taken 
  t_rx = 0.0 #corrected gps time
  dPR = [0]*NUM_CHANS #pseudo range error
  jac = [[0]*4 for _ in range(NUM_CHANS)]
  ma = [[0]*4 for _ in range(4)]
  mb = [[0]*4 for _ in range(4)]
  mc = [[0]*chans for _ in range(4)]
  md = [0]*4
  weight = [1]*NUM_CHANS

  p.x = p.y = p.z = p.time = t_pc = 0.0
  for i in range(chans):
      x_sv[i] = sv_l[i].x
      y_sv[i] = sv_l[i].y
      z_sv[i] = sv_l[i].z
      t_tx[i] = sv_l[i].time
      t_pc += sv_l[i].time
  #Approximate starting value for reciever clock
  t_pc = t_pc / chans + 75e-3 
  #iterate to user xyzt solution using "Taylor Series expansion"
  for j in range(MAX_ITER):
    t_rx = t_pc - p.time
    for i in range(chans):
      # Convert SV position to ECI coords(20.3.3.4.3.3.2)
      theta = (t_tx[i] - t_rx) * OMEGA_E

      x_sv_eci = x_sv[i]*math.cos(theta) - y_sv[i]*math.sin(theta)
      y_sv_eci = x_sv[i]*math.sin(theta) + y_sv[i]*math.cos(theta)
      z_sv_eci = z_sv[i]
      # Geometric range (20.3.3.4.3.4)
      gr = math.sqrt((p.x - x_sv_eci)**2 + (p.y - y_sv_eci)**2 + (p.z - z_sv_eci)**2)
      dPR[i] = SPEED_OF_LIGHT * (t_rx - t_tx[i]) - gr

      jac[i][0] = (p.x - x_sv_eci)/gr
      jac[i][1] = (p.y - y_sv_eci)/gr
      jac[i][2] = (p.z - z_sv_eci)/gr
      jac[i][3] = SPEED_OF_LIGHT

        #ma=transpose(H)*W*H
    ma = [[0 for c in range(4)] for r in range(4)]
    for r in range(4):
          for c in range(4):
              for i in range(chans):
                  ma[r][c] += jac[i][r] * weight[i] * jac[i][c]
    
    determinant = (
    ma[0][3]*ma[1][2]*ma[2][1]*ma[3][0] - ma[0][2]*ma[1][3]*ma[2][1]*ma[3][0] - ma[0][3]*ma[1][1]*ma[2][2]*ma[3][0] + ma[0][1]*ma[1][3]*ma[2][2]*ma[3][0] +ma[0][2]*ma[1][1]*ma[2][3]*ma[3][0] - ma[0][1]*ma[1][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[1][2]*ma[2][0]*ma[3][1] + ma[0][2]*ma[1][3]*ma[2][0]*ma[3][1] +ma[0][3]*ma[1][0]*ma[2][2]*ma[3][1] - ma[0][0]*ma[1][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[1][0]*ma[2][3]*ma[3][1] + ma[0][0]*ma[1][2]*ma[2][3]*ma[3][1] +ma[0][3]*ma[1][1]*ma[2][0]*ma[3][2] - ma[0][1]*ma[1][3]*ma[2][0]*ma[3][2] - ma[0][3]*ma[1][0]*ma[2][1]*ma[3][2] + ma[0][0]*ma[1][3]*ma[2][1]*ma[3][2] +ma[0][1]*ma[1][0]*ma[2][3]*ma[3][2] - ma[0][0]*ma[1][1]*ma[2][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[2][0]*ma[3][3] + ma[0][1]*ma[1][2]*ma[2][0]*ma[3][3] +ma[0][2]*ma[1][0]*ma[2][1]*ma[3][3] - ma[0][0]*ma[1][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[1][0]*ma[2][2]*ma[3][3] + ma[0][0]*ma[1][1]*ma[2][2]*ma[3][3])

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
                for i in range(4):
                    mc[r][c] +=mb[r][i]*jac[c][i]

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

  #return j

#************************************************************************
def LatLonAlt(x_n, y_n, z_n):
    iterations=100
    a = WGS84_A
    e2 = WGS84_E2

    p = math.sqrt(x_n*x_n + y_n*y_n)
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



#calling main function 
if __name__ == "__main__":
    main()

