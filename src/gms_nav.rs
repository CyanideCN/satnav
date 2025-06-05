use std::f64::consts::PI;

// Constants for coordinate conversions
const HALF_PI: f64 = PI / 2.0;

// Earth ellipsoid parameters (WGS84-like)
const SEMI_MAJOR_AXIS: f64 = 6_378_136.0; // meters
const FLATTENING: f64 = 1.0 / 298.257;
const FLATTENING_SQ: f64 = (1.0 - FLATTENING) * (1.0 - FLATTENING);
const ECCENTRICITY_SQ: f64 = 2.0 * FLATTENING - FLATTENING * FLATTENING;
/// 3D vector operations
/// Compute dot product of two 3D vectors
pub fn dot(vector_a: &[f64; 3], vector_b: &[f64; 3]) -> f64 {
    vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2]
}

/// Compute magnitude (length) of a 3D vector
pub fn magnitude(vector: &[f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

/// Normalize a 3D vector to unit length
pub fn normalize(vector: &[f64; 3]) -> [f64; 3] {
    let mag = magnitude(vector);
    if mag == 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [vector[0] / mag, vector[1] / mag, vector[2] / mag]
    }
}

/// Compute cross product of two 3D vectors
pub fn cross(vector_a: &[f64; 3], vector_b: &[f64; 3]) -> [f64; 3] {
    [
        vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1],
        vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2],
        vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0],
    ]
}

/// Compute normalized cross product of two 3D vectors
pub fn cross_normalized(vector_a: &[f64; 3], vector_b: &[f64; 3]) -> [f64; 3] {
    normalize(&cross(vector_a, vector_b))
}

/// Compute angle between two 3D vectors in radians
pub fn angle_between(vector_a: &[f64; 3], vector_b: &[f64; 3]) -> f64 {
    let denominator = magnitude(vector_a) * magnitude(vector_b);
    if denominator == 0.0 {
        0.0
    } else {
        (dot(vector_a, vector_b) / denominator).acos()
    }
}

/// Angular interpolation handling wrap-around
fn interp_angle(start_angle: f64, end_angle: f64, interp_factor: f64) -> f64 {
    // Both angles assumed to be in [0, 2π) range
    // Find shortest signed difference in (-π, π] range
    let mut angle_diff = end_angle - start_angle;

    // Handle wrap-around cases
    if angle_diff > PI {
        angle_diff -= 2.0 * PI; // wrap left
    }
    if angle_diff < -PI {
        angle_diff += 2.0 * PI; // wrap right
    }

    // Linear interpolation along shortest arc
    start_angle + angle_diff * interp_factor
}

/// Complete satellite state including position, orientation, and timing
#[derive(Debug)]
pub struct SatelliteState {
    pub position: [f64; 3],     // Satellite position in ECEF coordinates (meters)
    pub spin_axis: [f64; 3],    // Satellite spin axis unit vector
    pub sun_vec: [f64; 3],      // Sun direction unit vector
    pub beta_angle: f64,        // Beta angle for instrument alignment (radians)
    pub trans_matrix: [f64; 9], // 3x3 transformation matrix (row-major)
    pub sat_att_angle: f64,     // Satellite attitude angle (radians)
    pub sun_alpha: f64,         // Sun alpha angle (radians)
    pub sun_delta: f64,         // Sun delta angle (radians)
}

pub struct GMSNavigation {
    // Prediction tables for orbit and attitude
    orbit_pred_table: [[f64; 18]; 35],
    attitude_pred_table: [[f64; 33]; 10],

    // Instrument geometry parameters
    line_step_size: f64,         // Step size between scan lines (radians)
    pixel_sampling_rate: f64,    // Pixel sampling rate (radians per pixel)
    ref_line_center: f64,        // Reference line for center of scan (line number)
    ref_pixel_center: f64,       // Reference pixel for center of scan (pixel number)
    scan_start_time_offset: f64, // Time offset for scan start (seconds)
    spin_rate: f64,              // Satellite spin rate (Hz)
    sensor_lines_per_step: f64,  // Number of sensor lines per step

    // Sensor mounting and alignment
    sensor_misalignment: [f64; 3], // Sensor misalignment vector [x, y, z]
    sensor_mounting_mat: [[f64; 3]; 3], // 3x3 sensor mounting matrix
}

impl GMSNavigation {
    /// Create a new GMS Navigation instance
    pub fn new(
        attitude_pred_table: [[f64; 33]; 10],
        orbit_pred_table: [[f64; 18]; 35],
        line_step_size: f64,
        pixel_sampling_rate: f64,
        ref_line_center: f64,
        ref_pixel_center: f64,
        scan_start_time_offset: f64,
        spin_rate: f64,
        sensor_lines_per_step: f64,
        sensor_misalignment: [f64; 3],
        sensor_mounting_mat: [[f64; 3]; 3],
    ) -> Self {
        Self {
            orbit_pred_table,
            attitude_pred_table,
            line_step_size,
            pixel_sampling_rate,
            ref_line_center,
            ref_pixel_center,
            scan_start_time_offset,
            spin_rate,
            sensor_lines_per_step,
            sensor_misalignment,
            sensor_mounting_mat,
        }
    }

    /// Compute observation time relative to scan start for given line/pixel
    fn relative_observation_time(&self, line: f64, pixel: f64) -> f64 {
        let pixel_1_based = pixel + 1.0;

        let spin_frequency = 1440.0 * self.spin_rate; // lines per second
        let line_step = (line / self.sensor_lines_per_step).floor();
        let pixel_step = (self.pixel_sampling_rate * pixel_1_based) / (2.0 * PI);

        (line_step + pixel_step) / spin_frequency
    }

    /// Interpolate orbital and attitude data for a specific time interval
    fn interp_orbital_data(
        &self,
        table_index: usize,
        relative_time: f64,
    ) -> ([f64; 3], [f64; 9], f64, f64, f64) {
        let time_fraction = (relative_time - self.orbit_pred_table[0][table_index])
            / (self.orbit_pred_table[0][table_index + 1] - self.orbit_pred_table[0][table_index]);

        // Interpolate satellite position (ECEF coordinates)
        let sat_pos = [
            self.orbit_pred_table[8][table_index]
                + (self.orbit_pred_table[8][table_index + 1]
                    - self.orbit_pred_table[8][table_index])
                    * time_fraction,
            self.orbit_pred_table[9][table_index]
                + (self.orbit_pred_table[9][table_index + 1]
                    - self.orbit_pred_table[9][table_index])
                    * time_fraction,
            self.orbit_pred_table[10][table_index]
                + (self.orbit_pred_table[10][table_index + 1]
                    - self.orbit_pred_table[10][table_index])
                    * time_fraction,
        ];

        // Handle angular interpolation for satellite attitude angle
        let mut sat_att_angle = (self.orbit_pred_table[14][table_index]
            + (self.orbit_pred_table[14][table_index + 1]
                - self.orbit_pred_table[14][table_index])
                * time_fraction)
            .to_radians();

        if self.orbit_pred_table[14][table_index + 1] - self.orbit_pred_table[14][table_index] < 0.0
        {
            sat_att_angle = (self.orbit_pred_table[14][table_index]
                + (self.orbit_pred_table[14][table_index + 1]
                    - self.orbit_pred_table[14][table_index]
                    + 360.0)
                    * time_fraction)
                .to_radians();
        }

        // Handle angular interpolation for sun alpha angle
        let mut sun_alpha = (self.orbit_pred_table[17][table_index]
            + (self.orbit_pred_table[17][table_index + 1]
                - self.orbit_pred_table[17][table_index])
                * time_fraction)
            .to_radians();

        if self.orbit_pred_table[17][table_index + 1] - self.orbit_pred_table[17][table_index] > 0.0
        {
            sun_alpha = (self.orbit_pred_table[17][table_index]
                + (self.orbit_pred_table[17][table_index + 1]
                    - self.orbit_pred_table[17][table_index]
                    - 360.0)
                    * time_fraction)
                .to_radians();
        }

        // Interpolate sun delta angle
        let sun_delta = (self.orbit_pred_table[18][table_index]
            + (self.orbit_pred_table[18][table_index + 1]
                - self.orbit_pred_table[18][table_index])
                * time_fraction)
            .to_radians();

        // Extract 3x3 transformation matrix (stored in row-major order)
        let trans_matrix = [
            self.orbit_pred_table[19][table_index], // [0,0]
            self.orbit_pred_table[22][table_index], // [0,1]
            self.orbit_pred_table[25][table_index], // [0,2]
            self.orbit_pred_table[20][table_index], // [1,0]
            self.orbit_pred_table[23][table_index], // [1,1]
            self.orbit_pred_table[26][table_index], // [1,2]
            self.orbit_pred_table[21][table_index], // [2,0]
            self.orbit_pred_table[24][table_index], // [2,1]
            self.orbit_pred_table[27][table_index], // [2,2]
        ];

        (sat_pos, trans_matrix, sat_att_angle, sun_alpha, sun_delta)
    }

    /// Compute complete satellite state (position, spin axis, sun vector) for given time
    fn satellite_state(&self, relative_time: f64) -> SatelliteState {
        let mut sat_pos = [0.0; 3];
        let mut trans_mat = [0.0; 9];
        let mut sat_att_angle = 0.0;
        let mut sun_alpha = 0.0;
        let mut sun_delta = 0.0;

        // Find appropriate orbital data time interval and interpolate
        for i in 0..17 {
            if relative_time >= self.orbit_pred_table[0][i]
                && relative_time < self.orbit_pred_table[0][i + 1]
            {
                let (pos, trans_matrix, att_angle, s_alpha, s_delta) =
                    self.interp_orbital_data(i, relative_time);
                sat_pos = pos;
                trans_mat = trans_matrix;
                sat_att_angle = att_angle;
                sun_alpha = s_alpha;
                sun_delta = s_delta;
                break;
            }
        }

        // Interpolate attitude data from attitude prediction table
        let mut att_alpha = 0.0;
        let mut att_delta = 0.0;
        let mut beta_angle = 0.0;

        for i in 0..32 {
            if relative_time >= self.attitude_pred_table[0][i]
                && relative_time < self.attitude_pred_table[0][i + 1]
            {
                let time_fraction = (relative_time - self.attitude_pred_table[0][i])
                    / (self.attitude_pred_table[0][i + 1] - self.attitude_pred_table[0][i]);

                att_alpha = interp_angle(
                    self.attitude_pred_table[2][i],
                    self.attitude_pred_table[2][i + 1],
                    time_fraction,
                );
                att_delta = interp_angle(
                    self.attitude_pred_table[3][i],
                    self.attitude_pred_table[3][i + 1],
                    time_fraction,
                );
                beta_angle = interp_angle(
                    self.attitude_pred_table[4][i],
                    self.attitude_pred_table[4][i + 1],
                    time_fraction,
                );
                break;
            }
        }

        // Compute attitude vector in body frame
        let att_cos_delta = att_delta.cos();
        let att_vec = [
            att_delta.sin(),
            -att_cos_delta * att_alpha.sin(),
            att_cos_delta * att_alpha.cos(),
        ];

        // Transform attitude vector using transformation matrix
        let trans_att_vec = [
            trans_mat[0] * att_vec[0] + trans_mat[1] * att_vec[1] + trans_mat[2] * att_vec[2],
            trans_mat[3] * att_vec[0] + trans_mat[4] * att_vec[1] + trans_mat[5] * att_vec[2],
            trans_mat[6] * att_vec[0] + trans_mat[7] * att_vec[1] + trans_mat[8] * att_vec[2],
        ];

        // Compute final spin axis vector
        let (sin_att, cos_att) = sat_att_angle.sin_cos();
        let spin_axis = normalize(&[
            cos_att * trans_att_vec[0] + sin_att * trans_att_vec[1],
            -sin_att * trans_att_vec[0] + cos_att * trans_att_vec[1],
            trans_att_vec[2],
        ]);

        // Compute sun vector
        let (sun_delta_sin, sun_delta_cos) = sun_delta.sin_cos();
        let sun_vec = [
            sun_delta_cos * sun_alpha.cos(),
            sun_delta_cos * sun_alpha.sin(),
            sun_delta_sin,
        ];

        SatelliteState {
            position: sat_pos,
            spin_axis,
            sun_vec,
            beta_angle,
            trans_matrix: trans_mat,
            sat_att_angle,
            sun_alpha,
            sun_delta,
        }
    }

    /// Forward transform: Convert image coordinates (line, pixel) to geographic coordinates (lat, lon)
    pub fn image_to_geographic(&self, line: f64, pixel: f64) -> (f64, f64) {
        let observation_time =
            self.relative_observation_time(line, pixel) + self.scan_start_time_offset;
        let sat_state = self.satellite_state(observation_time);

        // Compute satellite body frame coordinate system
        let body_y_prime = cross_normalized(&sat_state.spin_axis, &sat_state.sun_vec);
        let body_x_prime = cross_normalized(&body_y_prime, &sat_state.spin_axis);

        // Apply beta angle rotation
        let (beta_sin, beta_cos) = sat_state.beta_angle.sin_cos();
        let body_x_rotated = [
            body_y_prime[0] * beta_sin + body_x_prime[0] * beta_cos,
            body_y_prime[1] * beta_sin + body_x_prime[1] * beta_cos,
            body_y_prime[2] * beta_sin + body_x_prime[2] * beta_cos,
        ];

        let body_x = normalize(&body_x_rotated);
        let body_y = cross_normalized(&sat_state.spin_axis, &body_x);

        // Compute instrument pointing angles
        let (line_angle_sin, line_angle_cos) =
            ((line - self.ref_line_center) * self.line_step_size).sin_cos();
        let (pixel_angle_sin, pixel_angle_cos) =
            ((pixel - self.ref_pixel_center) * self.pixel_sampling_rate).sin_cos();

        // Apply sensor mounting matrix transformation
        let sensor_dir_mat = [
            self.sensor_mounting_mat[0][0] * line_angle_cos
                + self.sensor_mounting_mat[0][2] * line_angle_sin,
            self.sensor_mounting_mat[1][0] * line_angle_cos
                + self.sensor_mounting_mat[1][2] * line_angle_sin,
            self.sensor_mounting_mat[2][0] * line_angle_cos
                + self.sensor_mounting_mat[2][2] * line_angle_sin,
        ];

        let sensor_dir_rotated = [
            pixel_angle_cos * sensor_dir_mat[0] - pixel_angle_sin * sensor_dir_mat[1],
            pixel_angle_sin * sensor_dir_mat[0] + pixel_angle_cos * sensor_dir_mat[1],
            sensor_dir_mat[2],
        ];

        // Transform to Earth-Centered Earth-Fixed (ECEF) coordinates
        let look_direction_ecef = [
            body_x[0] * sensor_dir_rotated[0]
                + body_y[0] * sensor_dir_rotated[1]
                + sat_state.spin_axis[0] * sensor_dir_rotated[2],
            body_x[1] * sensor_dir_rotated[0]
                + body_y[1] * sensor_dir_rotated[1]
                + sat_state.spin_axis[1] * sensor_dir_rotated[2],
            body_x[2] * sensor_dir_rotated[0]
                + body_y[2] * sensor_dir_rotated[1]
                + sat_state.spin_axis[2] * sensor_dir_rotated[2],
        ];

        let look_direction = normalize(&look_direction_ecef);

        // Compute intersection with Earth ellipsoid

        let ellipsoid_a = FLATTENING_SQ * (look_direction[0].powi(2) + look_direction[1].powi(2))
            + look_direction[2].powi(2);
        let ellipsoid_b = FLATTENING_SQ
            * (sat_state.position[0] * look_direction[0]
                + sat_state.position[1] * look_direction[1])
            + sat_state.position[2] * look_direction[2];
        let ellipsoid_c = FLATTENING_SQ
            * (sat_state.position[0].powi(2) + sat_state.position[1].powi(2)
                - SEMI_MAJOR_AXIS.powi(2))
            + sat_state.position[2].powi(2);

        let discriminant = ellipsoid_b.powi(2) - ellipsoid_a * ellipsoid_c;

        if discriminant <= 0.0 || ellipsoid_a == 0.0 {
            return (f64::NAN, f64::NAN);
        }

        let sqrt_discriminant = discriminant.sqrt();
        let intersect_dist_1 = (-ellipsoid_b + sqrt_discriminant) / ellipsoid_a;
        let intersect_dist_2 = (-ellipsoid_b - sqrt_discriminant) / ellipsoid_a;
        let intersect_dist = if intersect_dist_1.abs() <= intersect_dist_2.abs() {
            intersect_dist_1
        } else {
            intersect_dist_2
        };

        // Compute Earth intersection point
        let earth_point = [
            sat_state.position[0] + intersect_dist * look_direction[0],
            sat_state.position[1] + intersect_dist * look_direction[1],
            sat_state.position[2] + intersect_dist * look_direction[2],
        ];

        // Convert ECEF to geodetic coordinates
        let latitude_rad = (earth_point[2]
            / (FLATTENING_SQ * (earth_point[0].powi(2) + earth_point[1].powi(2)).sqrt()))
        .atan();
        let longitude_rad = if earth_point[0] != 0.0 {
            let mut lon = (earth_point[1] / earth_point[0]).atan();
            if earth_point[0] < 0.0 && earth_point[1] >= 0.0 {
                lon += PI;
            } else if earth_point[0] < 0.0 && earth_point[1] < 0.0 {
                lon -= PI;
            }
            lon
        } else if earth_point[1] > 0.0 {
            HALF_PI
        } else {
            -HALF_PI
        };

        let latitude = latitude_rad.to_degrees();
        let longitude = longitude_rad.to_degrees();
        (latitude, longitude)
    }

    /// Inverse transform: Convert geographic coordinates (lat, lon) to image coordinates (line, pixel)
    pub fn geographic_to_image_at_time(
        &self,
        longitude: f64,
        latitude: f64,
        relative_time: f64,
    ) -> (f64, f64) {
        let state = self.satellite_state(relative_time);

        // Convert geographic coordinates to Earth-Centered Earth-Fixed (ECEF)
        let latitude_rad = latitude.to_radians();
        let longitude_rad = longitude.to_radians();

        
        let prime_vertical_radius =
            SEMI_MAJOR_AXIS / (1.0 - ECCENTRICITY_SQ * latitude_rad.sin().powi(2)).sqrt();

        let earth_point = [
            prime_vertical_radius * latitude_rad.cos() * longitude_rad.cos(),
            prime_vertical_radius * latitude_rad.cos() * longitude_rad.sin(),
            prime_vertical_radius * (1.0 - ECCENTRICITY_SQ) * latitude_rad.sin(),
        ];

        // Compute look vector from satellite to Earth point
        let look_vector_raw = [
            earth_point[0] - state.position[0],
            earth_point[1] - state.position[1],
            earth_point[2] - state.position[2],
        ];
        let look_vector = normalize(&look_vector_raw);

        // Compute satellite body frame coordinate system
        let body_y_prime = cross_normalized(&state.spin_axis, &state.sun_vec);
        let body_x_prime = cross_normalized(&body_y_prime, &state.spin_axis);

        let (beta_sin, beta_cos) = state.beta_angle.sin_cos();
        let body_x_rotated = [
            body_y_prime[0] * beta_sin + body_x_prime[0] * beta_cos,
            body_y_prime[1] * beta_sin + body_x_prime[1] * beta_cos,
            body_y_prime[2] * beta_sin + body_x_prime[2] * beta_cos,
        ];
        let body_x = normalize(&body_x_rotated);
        let body_y = cross_normalized(&state.spin_axis, &body_x);

        // Compute working angles for instrument pointing
        let pixel_angle = {
            let cross_spin_look = cross(&state.spin_axis, &look_vector);
            let cross_y_cross = cross(&body_y, &cross_spin_look);
            let mut angle = angle_between(&body_y, &cross_spin_look);
            if dot(&state.spin_axis, &cross_y_cross) < 0.0 {
                angle = -angle;
            }
            angle
        };
        let line_angle = angle_between(&state.spin_axis, &look_vector);

        // Map to detector coordinates
        let line = (HALF_PI - line_angle) / self.line_step_size + self.ref_line_center
            - self.sensor_misalignment[1] / self.line_step_size;
        let pixel = pixel_angle / self.pixel_sampling_rate
            + self.ref_pixel_center
            + self.sensor_misalignment[2] / self.pixel_sampling_rate
            - (HALF_PI - line_angle) * self.sensor_misalignment[0].tan() / self.pixel_sampling_rate;

        (line, pixel)
    }

    pub fn geographic_to_image(&self, longitude: f64, latitude: f64) -> (f64, f64) {
        let (line_approx, pixel_approx) =
            self.geographic_to_image_at_time(longitude, latitude, self.scan_start_time_offset);
        // Use image coords computed at scan start time to compute a more accurate relative time
        let relative_time =
            self.relative_observation_time(line_approx, pixel_approx) + self.scan_start_time_offset;
        let (line, pixel) = self.geographic_to_image_at_time(longitude, latitude, relative_time);
        (line, pixel)
    }
}
