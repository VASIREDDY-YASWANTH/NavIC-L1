
import math

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




class Location:
    def __init__(self, X , Y , Z, T):
        self.x = X
        self.y = Y
        self.z = Z
        self.time =T

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


S_V=[]
S_V.append(Space_vehicle(1, 0.627712271718553, 5153.52623939514, 6.60E-06, -3.65E-06, 1.30E-08, 1.30E-08, 254.875, -73.15625, 4.56E-09, 0.00440125516615808, 0.948704618521089, 0.543838107760104, -5.75E-11, -8.17E-09, 0.96682775824546, -1.91E-08, 468000, 0, -3.64E-12, -0.000181447714567184, 0, 468000,    873300956  , 0  , 77786, 224, 0))
S_V.append(Space_vehicle(2, 2.89655783543179, 5153.61238670349, 7.14E-06, 6.39E-07, -2.22E-07, -5.59E-08, 224.5625, 16, 5.20E-09, 0.0169517478207126, -2.43433635302316, 2.65227448464399, 1.57E-11, -8.06E-09, 0.932250669788079, -1.16E-08, 468000 , 0, -2.73E-12, -0.000101803801953793, 0, 468000,    4099961621   , 0,   77786  , 218, 0)) 
S_V.append(Space_vehicle(3, 0.927649448643036, 5153.74806022644, 6.28E-06, -3.84E-06, 5.59E-09, 1.12E-08, 257.46875, -74.40625, 4.59E-09, 0.000642508151941001, 1.94139717071046, 0.531753256034665, -2.75E-11, -8.13E-09, 0.961676863187149, 4.66E-09, 467952, 0, -2.39E-12, -5.12E-05, 0, 467952,   1716325164  , 0,   77786 , 218, 0)) 
S_V.append(Space_vehicle(4, 3.02458803990085, 5153.63453292847, 9.39E-06, -1.34E-06, 4.84E-08, 5.22E-08, 196.9375, -25.53125, 4.49E-09, 0.00251956866122782, -1.25409232099612, 1.62825900432321, 1.35E-10,  -7.87E-09, 0.960616162730431, -8.38E-09, 468000, 0, 2.84E-12, 1.63E-04, 0, 468000,   100389401  , 0,  77786 , 225,  0  )) 
S_V.append(Space_vehicle(5, 0.796297202986197, 5155.34639358521, 6.45E-06, -3.75E-06, -6.89E-08, -1.60E-07, 256.9375, -70.65625, 4.81E-09, 0.0130847115069628, 1.55579952931809, 0.503145066445881, -5.89E-11, -8.41E-09, 0.958335540739085, -7.92E-09, 468000, 0, 3.41E-12, 2.96E-04, 0, 468000,  3531128408 
  , 0,    77786   , 226,  0 )) 
S_V.append(Space_vehicle(6, -2.60055087716526, 5153.80060768127, 4.23E-06, 3.95E-07, 1.30E-07, -1.49E-08, 305.5625, 3.8125, 4.28E-09, 0.00755723693873733, -1.04896671410838, -0.494186142962186, 1.69E-10, -7.95E-09, 0.978481871297476, -1.30E-08, 467984, 0, 3.75E-12, 7.98E-06, 0, 467984,  2182195956 
 , 0,   77786  , 223, 0))  
L=[]
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))
L.append(Location(0.0 , 0.0 , 0.0 , 0.0))

def sv_calc_corrected_time(sv):
  ##  print("               cala_corrected_time   ")
    #global Space_vehicle
    #sv = Space_vehicle(i)
    #global s_v
    #sv = s_v[i]

    delta_t = sv.time_raw - sv.reference_time
    #print("time raw= \n",sv.time_raw)
    #print("reference time= \n",sv.reference_time)
    if(delta_t >  302400):  
        delta_t -= 604800
    if(delta_t < -302400):  
        delta_t += 604800

    ek = orbit_ecc_anom(sv, sv.time_raw)

    #print(" ek= \n",ek)
    #print(" e= \n",sv.e)
    #print("sqrt sqrtA= \n",sv.sqrt_A)
    delta_tr = -4.442807633e-10 * sv.e * sv.sqrt_A * math.sin(ek)

    time_correction = sv.correction_f0 + (sv.correction_f1 * delta_t) + (sv.correction_f2 * delta_t * delta_t) + delta_tr - sv.group_delay
    #print("time_correction= \n",time_correction)
    sv.pos_t = sv.time_raw - time_correction
    ##print("pos_t ",sv.pos_t,"\n")




def orbit_ecc_anom(sv, t):
    ##print("               orbit_ecc_anom   ")
    iterations  = 200
    delta       = 1e-10
    #print("sqrt a= \n",sv.sqrt_A)
    semi_major_axis = sv.sqrt_A * sv.sqrt_A
    computed_mean_motion = math.sqrt(mu / pow(semi_major_axis,3.0))
    #print("time of epi= \n",sv.time_of_ephemeris)
    time_from_ephemeris  = t - sv.time_of_ephemeris
    if(time_from_ephemeris >  302400):
          time_from_ephemeris -= 604800
    if(time_from_ephemeris < -302400):
          time_from_ephemeris += 604800
    #print("delta n= \n",sv.delta_n)
    correct_mean_motion  = computed_mean_motion + sv.delta_n

    #print("correct mean motion = \n",correct_mean_motion)
    # Add on how much we have moved through the orbit since ephemeris
    print("mean mot at epi= \n",sv.mean_motion_at_ephemeris)
    mean_anomaly = sv.mean_motion_at_ephemeris + correct_mean_motion * time_from_ephemeris
    #print("mean anomoly= \n",mean_anomaly)
    estimate = ( mean_anomaly if sv.e < 0.8 else math.pi)
    correction = estimate - (mean_anomaly + sv.e*math.sin(mean_anomaly))

    # Solve iteratively
    while ((abs(correction) > delta) and iterations > 0):
        last = estimate
        estimate  = mean_anomaly  + sv.e * math.sin(estimate)
        correction = estimate - last
        iterations -= 1
    #print("estimate= \n",estimate)
    return estimate

def sv_calc_location(sv, l: Location):
    l.time = sv.pos_t
    ##print("l->time= \n",l.time)
    orbit_calc_position(sv, l)
    return 1

def orbit_calc_position(sv, l: Location):
    time_from_ephemeris = sv.pos_t - sv.time_of_ephemeris
    if time_from_ephemeris > 302400:
        time_from_ephemeris -= 604800
    if time_from_ephemeris < -302400:
        time_from_ephemeris += 604800

    semi_major_axis = sv.sqrt_A * sv.sqrt_A
    ek = orbit_ecc_anom(sv, sv.pos_t)
    ##print("ek in orb pos cal= \n",ek)
    # Now calculate the first approximation of the latitude
    true_anomaly = math.atan2(math.sqrt(1 - sv.e * sv.e) * math.sin(ek), math.cos(ek) - sv.e)
    ##print("true_anomaly= \n",true_anomaly)
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
    print("l->x= \n",l.x)
    print("l->y= \n",l.y)
    print("l->z= \n",l.z)
    print("l->time= \n",l.time)
    return 1



def A_solve(chans, sv_l , p ):  
  print("            A solve             ") 
  print("chans= \n",chans)
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
      x_sv[i] = sv_l[i].x
      y_sv[i] = sv_l[i].y
      z_sv[i] = sv_l[i].z
      t_tx[i] = sv_l[i].time
      t_pc += sv_l[i].time
      print("i,sv_l \n",i, sv_l[i].x, sv_l[i].y, sv_l[i].z, sv_l[i].time )
  t_pc = t_pc / chans + 75e-3 

  print("t_pc = \n",t_pc)
  for j in range(MAX_ITER):
  #for j in range(5):
    t_rx = t_pc - p.time
    for i in range(chans):
      # Convert SV position to ECI coords
      theta = (t_tx[i] - t_rx) * OMEGA_E

      x_sv_eci = x_sv[i]*math.cos(theta) - y_sv[i]*math.sin(theta)
      y_sv_eci = x_sv[i]*math.sin(theta) + y_sv[i]*math.cos(theta)
      z_sv_eci = z_sv[i]
      #print("i,eci \n",i, x_sv_eci, y_sv_eci, z_sv_eci )
      # Geometric range 
      gr = math.sqrt((p.x - x_sv_eci)**2 + (p.y - y_sv_eci)**2 + (p.z - z_sv_eci)**2)

      #print("gr= \n",gr)
      dPR[i] = SPEED_OF_LIGHT * (t_rx - t_tx[i]) - gr

      jac[i][0] = (p.x - x_sv_eci)/gr
      jac[i][1] = (p.y - y_sv_eci)/gr
      jac[i][2] = (p.z - z_sv_eci)/gr
      jac[i][3] = SPEED_OF_LIGHT

      #print("i,jac \n",i, jac[i][0],jac[i][1],jac[i][2],jac[i][3] )
        #ma=transpose(H)*W*H
    ma = [[0 for c in range(4)] for r in range(4)]
    for r in range(4):
          for c in range(4):
              for i in range(chans):
                  ma[r][c] += jac[i][r] * weight[i] * jac[i][c]
    '''  for r in range(4):
          for c in range(4):
              print(ma[r][c], "\t",end="" )
          print("\n")
    '''
    determinant = (
    ma[0][3]*ma[1][2]*ma[2][1]*ma[3][0] - ma[0][2]*ma[1][3]*ma[2][1]*ma[3][0] - ma[0][3]*ma[1][1]*ma[2][2]*ma[3][0] + ma[0][1]*ma[1][3]*ma[2][2]*ma[3][0] +ma[0][2]*ma[1][1]*ma[2][3]*ma[3][0] - ma[0][1]*ma[1][2]*ma[2][3]*ma[3][0] - ma[0][3]*ma[1][2]*ma[2][0]*ma[3][1] + ma[0][2]*ma[1][3]*ma[2][0]*ma[3][1] +ma[0][3]*ma[1][0]*ma[2][2]*ma[3][1] - ma[0][0]*ma[1][3]*ma[2][2]*ma[3][1] - ma[0][2]*ma[1][0]*ma[2][3]*ma[3][1] + ma[0][0]*ma[1][2]*ma[2][3]*ma[3][1] +ma[0][3]*ma[1][1]*ma[2][0]*ma[3][2] - ma[0][1]*ma[1][3]*ma[2][0]*ma[3][2] - ma[0][3]*ma[1][0]*ma[2][1]*ma[3][2] + ma[0][0]*ma[1][3]*ma[2][1]*ma[3][2] +ma[0][1]*ma[1][0]*ma[2][3]*ma[3][2] - ma[0][0]*ma[1][1]*ma[2][3]*ma[3][2] - ma[0][2]*ma[1][1]*ma[2][0]*ma[3][3] + ma[0][1]*ma[1][2]*ma[2][0]*ma[3][3] +ma[0][2]*ma[1][0]*ma[2][1]*ma[3][3] - ma[0][0]*ma[1][2]*ma[2][1]*ma[3][3] - ma[0][1]*ma[1][0]*ma[2][2]*ma[3][3] + ma[0][0]*ma[1][1]*ma[2][2]*ma[3][3])
    print("determinant\n",determinant)
    determinant = 5542836919226368
    print("determinant\n",determinant,"\n")

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

    #if err_mag < 1.0:
      #break

    p.x += dx
    p.y += dy
    p.z += dz
    p.time += dt

  #return j


def LatLonAlt(x_n, y_n, z_n):
    print()
    #x_n=3851769.812880
    #y_n=-78310.313391
    #z_n=5066279.257806
    iterations=100
    a = WGS84_A
    e2 = WGS84_E2

    p = math.sqrt(x_n*x_n + y_n*y_n)
    print(" p \n",p)
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
    predicted_location = Location(0,0,0,0)
    sv_location = Location(0,0,0,0)
    valid_locations = 0
    sv_ids = [1, 2, 3, 4 , 5,6  ]
    size = len(sv_ids)
    #size = 1
    
    for i in range(size):
        print("********************Main Code********************\n")
        sv = S_V[i]  # Initialize SpaceVehicle object, replace with actual values
        sv_location=L[i]
        print("sv_id=\n",i)
        #sv.lock_code_nco=  3425404228
        #sv.lock_code_nco=0xD10E5D7C
        print("lock_code_nco",sv.lock_code_nco,"\n")
        phase_gold_code = (sv.lock_code_nco >> 1) & 0x7FFFFFFF
        print("phase_gold_code =\n",phase_gold_code)
        #sv.nav_subframe_week=  77789
        #sv.lock_ms_frame= 66
        print("nav subf week =\n",sv.nav_subframe_week)
        print("lock_ms_frame =\n",sv.lock_ms_frame)
        sv.time_raw = sv.nav_subframe_week * 6000 + sv.lock_ms_frame + (float(phase_gold_code) / (1023 * (1 << 21)))
        sv.time_raw /= 1000.0
        print("time_raw in main code",sv.time_raw,"\n")
        sv_calc_corrected_time(sv)
        sv_calc_location(sv, sv_location)
        if -40000000 < sv_location.x < 40000000 and -40000000 < sv_location.y < 40000000 and -40000000 < sv_location.z < 40000000:
            valid_locations += 1

        if valid_locations < 4:
            sv_location
            #print(valid_locations,sv_location.x,sv_location.y,sv_location.z,sv_location.time)
            print(predicted_location.x,predicted_location.y,predicted_location.z,predicted_location.time)
        A_solve(valid_locations,L , predicted_location)
        lat,lon,alt = LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z)
        print(f"Lat/Lon/Alt : {lat * 180 / PI:20.6f}, {lon * 180 / PI:20.6f}, {alt:20.1f}")


if __name__ == "__main__":
    main()




















'''        
        sv_calc_corrected_time(i)
        sv_calc_location(i, sv_location+valid_locations)
        if -40000000 < sv_location.x < 40000000 and -40000000 < sv_location.y < 40000000 and -40000000 < sv_location.z < 40000000:
            valid_locations += 1

    if valid_locations < 4:
        sv_location.clear()

    A_solve(valid_locations, sv_location, predicted_location)
    LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z,lat,lon,alt)

    print(f"Lat/Lon/Alt : {lat * 180 / PI:20.6f}, {lon * 180 / PI:20.6f}, {alt:20.1f}")




def sv_calc_corrected_time(i):
    #global Space_vehicle
    #sv = Space_vehicle(i)
    global s_v
    sv = s_v[i]

    delta_t = sv.time_raw - sv.reference_time
    print(sv.time_raw)
    if(delta_t >  302400):  
        delta_t -= 604800
    if(delta_t < -302400):  
        delta_t += 604800

    ek = orbit_ecc_anom(sv, sv.time_raw)

    delta_tr = -4.442807633e-10 * sv.e * sv.sqrt_A * math.sin(ek)

    time_correction = sv.correction_f0 + (sv.correction_f1 * delta_t) + (sv.correction_f2 * delta_t * delta_t) + delta_tr - sv.group_delay
    sv.pos_t = sv.time_raw - time_correction 




if __name__ == "__main__":
    main()



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

'''

