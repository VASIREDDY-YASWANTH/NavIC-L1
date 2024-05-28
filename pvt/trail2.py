from typing import List, Dict

class SpaceVehicle:
    def __init__(self, ids):
        self.ids = ids
        self.mean_motion_at_ephemeris = None  # Initialize all variables
        self.sqrt_A = None
        self.Cus = None
        self.Cuc = None
        self.Cis = None
        self.Cic = None
        self.Crc = None
        self.Crs = None
        self.delta_n = None
        self.e = None
        self.w = None
        self.omega_0 = None
        self.idot = None
        self.omega_dot = None
        self.inclination_at_ephemeris = None
        self.group_delay = None
        self.reference_time = None
        self.correction_f2 = None
        self.correction_f1 = None
        self.correction_f0 = None
        self.pos_t = None
        self.time_correction = None
        self.time_of_ephemeris = None
        self.lock_code_nco = None
        self.time_raw = None
        self.nav_subframe_week = None
        self.lock_ms_frame = None
        self.subframe: List[List[int]] = None
        self.subframe_number = None
        self.synced = None

# Create a list of SpaceVehicle objects
space_vehicles: List[SpaceVehicle] = []
for ids in range(1, 33):  # Loop through 32 IDs
    space_vehicles.append(SpaceVehicle(ids))

class Location:
    def __init__(self, x, y, z, time):
        self.x = x
        self.y = y
        self.z = z
        self.time = time


from math import pi, sqrt

def orbit_ecc_anom(sv, t):
    """
    Calculates the eccentric anomaly of a satellite using iterative method.

    Args:
        sv: A SpaceVehicle object containing the necessary attributes.
        t: The time in seconds since epoch.

    Returns:
        The eccentric anomaly in radians.
    """

    iterations = 200
    delta = 1e-10
    estimate = None
    correction = None

    # Calculate semi-major axis and mean motion
    semi_major_axis = sv.sqrt_A ** 2
    computed_mean_motion = sqrt(mu / semi_major_axis ** 3)

    # Calculate time from ephemeris and adjust for week rollover
    time_from_ephemeris = t - sv.time_of_ephemeris
    if time_from_ephemeris > 302400:
        time_from_ephemeris -= 604800
    elif time_from_ephemeris < -302400:
        time_from_ephemeris += 604800

    # Calculate corrected mean motion
    correct_mean_motion = computed_mean_motion + sv.delta_n

    # Calculate mean anomaly
    mean_anomaly = sv.mean_motion_at_ephemeris + correct_mean_motion * time_from_ephemeris

    # Initial estimate
    if sv.e < 0.8:
        estimate = mean_anomaly
    else:
        estimate = pi

    # Iteratively solve for eccentric anomaly
    for _ in range(iterations):
        last_estimate = estimate
        estimate = mean_anomaly + sv.e * sin(estimate)
        correction = estimate - last_estimate
        if abs(correction) < delta:
            break

    return estimate




def sv_calc_corrected_time(sv):
    """
    Calculates the corrected time for a given SpaceVehicle object.

    Args:
        sv: A SpaceVehicle object containing the necessary attributes.

    Returns:
        The corrected time in seconds.
    """

    delta_t = sv.time_raw - sv.reference_time

    # Wrap delta_t to account for week rollover
    if delta_t > 302400:
        delta_t -= 604800
    elif delta_t < -302400:
        delta_t += 604800

    # Calculate eccentric anomaly (replace with actual implementation)
    ek = orbit_ecc_anom(sv, sv.time_raw)

    # Calculate relativistic term (replace with actual implementation)
    delta_tr = -4.442807633e-10 * sv.e * sv.sqrt_A * math.sin(ek)

    # Calculate time correction
    time_correction = (
        sv.correction_f0
        + sv.correction_f1 * delta_t
        + sv.correction_f2 * delta_t ** 2
        + delta_tr
        - sv.group_delay
    )

    # Calculate corrected position time
    sv.pos_t = sv.time_raw - time_correction

    return sv.pos_t


def sv_calc_location(sv_id, location):
    """
    Calculates the position of a satellite given its ID and updates a Location object.

    Args:
        sv_id: The ID of the satellite in the space_vehicles list.
        location: A Location object to be updated with the calculated position.

    Returns:
        None
    """

    sv = space_vehicles[sv_id]
    location.time = sv.pos_t

    # Call the corresponding function to calculate position
    orbit_calc_position(sv, location)


from math import atan2, cos, sin, sqrt

def orbit_calc_position(sv, location):
    """
    Calculates the position of a satellite given its state and updates a Location object.

    Args:
        sv: A SpaceVehicle object containing the necessary attributes.
        location: A Location object to be updated with the calculated position.

    Returns:
        None
    """

    # Time from ephemeris (adjust for week rollover)
    time_from_ephemeris = sv.pos_t - sv.time_of_ephemeris
    if time_from_ephemeris > 302400:
        time_from_ephemeris -= 604800
    elif time_from_ephemeris < -302400:
        time_from_ephemeris += 604800

    # Semi-major axis
    semi_major_axis = sv.sqrt_A * sv.sqrt_A

    # Eccentric anomaly (replace with actual implementation)
    ek = orbit_ecc_anom(sv, sv.pos_t)

    # True anomaly
    true_anomaly = atan2(sqrt(1 - sv.e * sv.e) * sin(ek), cos(ek) - sv.e)

    # Argument of latitude
    argument_of_latitude = true_anomaly + sv.w

    # Second harmonic perturbations
    argument_of_latitude_correction = (
        sv.Cus * sin(2 * argument_of_latitude) + sv.Cuc * cos(2 * argument_of_latitude)
    )
    radius_correction = (
        sv.Crc * cos(2 * argument_of_latitude) + sv.Crs * sin(2 * argument_of_latitude)
    )
    correction_of_inclination = (
        sv.Cic * cos(2 * argument_of_latitude) + sv.Cis * sin(2 * argument_of_latitude)
    )

    # Corrected values
    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
    corrected_radius = semi_major_axis * (1 - sv.e * cos(ek)) + radius_correction
    corrected_inclination = (
        sv.inclination_at_ephemeris + correction_of_inclination + sv.idot * time_from_ephemeris
    )

    # Position in orbital plane
    pos_in_orbial_plane_x = corrected_radius * cos(corrected_argument_of_latitude)
    pos_in_orbial_plane_y = corrected_radius * sin(corrected_argument_of_latitude)

    # Corrected angle of ascending node
    corrected_angle_of_ascending_node = (
        sv.omega_0
        + (sv.omega_dot - OMEGA_E) * time_from_ephemeris
        - OMEGA_E * sv.time_of_ephemeris
    )

    # Earth-centered, Earth-fixed coordinates
    location.x = (
        pos_in_orbial_plane_x * cos(corrected_angle_of_ascending_node)
        - pos_in_orbial_plane_y * cos(corrected_inclination) * sin(corrected_angle_of_ascending_node)
    )
    location.y = (
        pos_in_orbial_plane_x * sin(corrected_angle_of_ascending_node)
        + pos_in_orbial_plane_y * cos(corrected_inclination) * cos(corrected_angle_of_ascending_node)
    )
    location.z = pos_in_orbial_plane_y * sin(corrected_inclination)
    location.time = sv.pos_t

import math

NUM_CHANS = 12
MAX_ITER = 10
SPEED_OF_LIGHT = 299792458.0  # m/s
OMEGA_E = 7.292115e-5  # Earth rotation rate (rad/s)

def A_solve(chans, sv_l, p):
    """Iterates to user xyzt solution using Taylor Series expansion."""

    t_tx = [0.0] * chans
    x_sv = [0.0] * chans
    y_sv = [0.0] * chans
    z_sv = [0.0] * chans
    weight = [1.0] * chans

    p.x = p.y = p.z = p.time = t_pc = 0.0

    for i in range(chans):
        if i < NUM_CHANS:
            x_sv[i] = sv_l[i].x
            y_sv[i] = sv_l[i].y
            z_sv[i] = sv_l[i].z
            t_tx[i] = sv_l[i].time
            t_pc += sv_l[i].time

    t_pc /= chans + 75e-3

    for j in range(MAX_ITER):
        t_rx = t_pc - p.time

        jac = [[0.0] * 4 for _ in range(chans)]
        dPR = [0.0] * chans

        for i in range(chans):
            theta = (t_tx[i] - t_rx) * OMEGA_E

            x_sv_eci = x_sv[i] * math.cos(theta) - y_sv[i] * math.sin(theta)
            y_sv_eci = x_sv[i] * math.sin(theta) + y_sv[i] * math.cos(theta)
            z_sv_eci = z_sv[i]

            gr = math.sqrt(
                (p.x - x_sv_eci)**2 + (p.y - y_sv_eci)**2 + (p.z - z_sv_eci)**2
            )

            dPR[i] = SPEED_OF_LIGHT * (t_rx - t_tx[i]) - gr

            jac[i][0] = (p.x - x_sv_eci) / gr
            jac[i][1] = (p.y - y_sv_eci) / gr
            jac[i][2] = (p.z - z_sv_eci) / gr
            jac[i][3] = SPEED_OF_LIGHT

      
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


        # Matrix calculations (ma, mb, mc, md) are equivalent to the C code
        # ... (omitted for brevity)

        dx = md[0]
        dy = md[1]
        dz = md[2]
        dt = md[3]

        err_mag = math.sqrt(dx**2 + dy**2 + dz**2)

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
    lat, lon, alt = None, None, None  # Initialize variables
    predicted_location = Location(None, None, None, None)
    sv_location = [None] * 5  # Placeholder for calculated locations
    valid_locations = 0
    sv_ids = [1, 2, 3, 4, 5]  # Assuming satellite IDs

    for sv_id in sv_ids:
        sv = space_vehicles[sv_id - 1]

        # Calculate corrected time (replace with actual implementation)
        sv.time_raw = sv_calc_corrected_time(sv)  # Replace with actual function

        # Calculate satellite location (replace with actual implementation)
        sv_calc_location(sv_id - 1, sv_location)  # Replace with actual function

        if -40000000 < sv_location[valid_locations].x < 40000000 and \
           -40000000 < sv_location[valid_locations].y < 40000000 and \
           -40000000 < sv_location[valid_locations].z < 40000000:
            valid_locations += 1

    if valid_locations < 4:
        # Handle insufficient valid locations (e.g., raise error or return early)
        return

    # Solve for position using A_solve (replace with actual implementation)
    A_solve(valid_locations, sv_location, predicted_location)

    # Convert Cartesian to latitude/longitude/altitude (replace with actual implementation)
    lat, lon, alt = LatLonAlt(predicted_location.x, predicted_location.y, predicted_location.z)

    print(f"Lat/Lon/Alt: {lat:.6f}, {lon:.6f}, {alt:.1f}")

if __name__ == "__main__":
    main()

