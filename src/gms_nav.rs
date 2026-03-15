use std::f64::consts::PI;

use crate::math3d::{
    add_scaled, angle_between, change_basis, cross, cross_normalized, dot, mat3_mul_vec3,
    normalize, subtract, Mat3, Vec3,
};

// Constants for coordinate conversions
const HALF_PI: f64 = PI / 2.0;

// Earth ellipsoid parameters (WGS84-like)
const SEMI_MAJOR_AXIS: f64 = 6_378_136.0; // meters
const FLATTENING: f64 = 1.0 / 298.257;
const FLATTENING_SQ: f64 = (1.0 - FLATTENING) * (1.0 - FLATTENING);
const ECCENTRICITY_SQ: f64 = 2.0 * FLATTENING - FLATTENING * FLATTENING;

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
    pub position: Vec3,     // Satellite position in ECEF coordinates (meters)
    pub spin_axis: Vec3,    // Satellite spin axis unit vector
    pub sun_vec: Vec3,      // Sun direction unit vector
    pub beta_angle: f64,    // Beta angle for instrument alignment (radians)
    pub conv_matrix: Mat3,  // 3x3 conversion matrix
    pub sat_att_angle: f64, // Satellite attitude angle (radians)
    pub sun_alpha: f64,     // Sun alpha angle (radians)
    pub sun_delta: f64,     // Sun delta angle (radians)
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
    misalignment_mat: [[f64; 3]; 3], // 3x3 misalignment matrix
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
        misalignment_mat: [[f64; 3]; 3],
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
            misalignment_mat,
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
    ) -> (Vec3, Mat3, f64, f64, f64) {
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

        // Extract 3x3 conversion matrix (stored in row-major order)
        let conv_matrix = [
            [
                self.orbit_pred_table[19][table_index],
                self.orbit_pred_table[22][table_index],
                self.orbit_pred_table[25][table_index],
            ],
            [
                self.orbit_pred_table[20][table_index],
                self.orbit_pred_table[23][table_index],
                self.orbit_pred_table[26][table_index],
            ],
            [
                self.orbit_pred_table[21][table_index],
                self.orbit_pred_table[24][table_index],
                self.orbit_pred_table[27][table_index],
            ],
        ];

        (sat_pos, conv_matrix, sat_att_angle, sun_alpha, sun_delta)
    }

    fn satellite_body_axes(
        &self,
        spin_axis: &Vec3,
        sun_vec: &Vec3,
        beta_angle: f64,
    ) -> (Vec3, Vec3) {
        let sat_y_prime = cross_normalized(spin_axis, sun_vec);
        let sat_x_prime = cross_normalized(&sat_y_prime, spin_axis);

        let (beta_sin, beta_cos) = beta_angle.sin_cos();
        let sat_x = normalize(&add_scaled(&sat_y_prime, beta_sin, &sat_x_prime, beta_cos));
        let sat_y = cross_normalized(spin_axis, &sat_x);

        (sat_x, sat_y)
    }

    /// Compute complete satellite state (position, spin axis, sun vector) for given time
    fn satellite_state(&self, relative_time: f64) -> SatelliteState {
        let mut sat_pos = [0.0; 3];
        let mut conv_mat = [[0.0; 3]; 3];
        let mut sat_att_angle = 0.0;
        let mut sun_alpha = 0.0;
        let mut sun_delta = 0.0;

        // Find appropriate orbital data time interval and interpolate
        for i in 0..17 {
            if relative_time >= self.orbit_pred_table[0][i]
                && relative_time < self.orbit_pred_table[0][i + 1]
            {
                let (pos, conv_matrix, att_angle, s_alpha, s_delta) =
                    self.interp_orbital_data(i, relative_time);
                sat_pos = pos;
                conv_mat = conv_matrix;
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

        // Transform attitude vector using conversion matrix
        let att_vec_conv = mat3_mul_vec3(&conv_mat, &att_vec);

        // Compute final spin axis vector
        let (sin_att, cos_att) = sat_att_angle.sin_cos();
        let spin_axis = normalize(&[
            cos_att * att_vec_conv[0] + sin_att * att_vec_conv[1],
            -sin_att * att_vec_conv[0] + cos_att * att_vec_conv[1],
            att_vec_conv[2],
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
            conv_matrix: conv_mat,
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

        let (sat_x, sat_y) = self.satellite_body_axes(
            &sat_state.spin_axis,
            &sat_state.sun_vec,
            sat_state.beta_angle,
        );

        // Compute instrument pointing angles
        let (line_angle_sin, line_angle_cos) =
            ((line - self.ref_line_center) * self.line_step_size).sin_cos();
        let (pixel_angle_sin, pixel_angle_cos) =
            ((pixel - self.ref_pixel_center) * self.pixel_sampling_rate).sin_cos();

        // Apply misalignment matrix
        let view_vec_sat = mat3_mul_vec3(
            &self.misalignment_mat,
            &[line_angle_cos, 0.0, line_angle_sin],
        );

        let view_vec_sat_rotated = [
            pixel_angle_cos * view_vec_sat[0] - pixel_angle_sin * view_vec_sat[1],
            pixel_angle_sin * view_vec_sat[0] + pixel_angle_cos * view_vec_sat[1],
            view_vec_sat[2],
        ];

        // Transform to Earth-Centered Earth-Fixed (ECEF) coordinates
        let view_vec_ecef =
            change_basis(&sat_x, &sat_y, &sat_state.spin_axis, &view_vec_sat_rotated);

        let view_vec = normalize(&view_vec_ecef);

        // Compute intersection with Earth ellipsoid

        let ellipsoid_a =
            FLATTENING_SQ * (view_vec[0].powi(2) + view_vec[1].powi(2)) + view_vec[2].powi(2);
        let ellipsoid_b = FLATTENING_SQ
            * (sat_state.position[0] * view_vec[0] + sat_state.position[1] * view_vec[1])
            + sat_state.position[2] * view_vec[2];
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
            sat_state.position[0] + intersect_dist * view_vec[0],
            sat_state.position[1] + intersect_dist * view_vec[1],
            sat_state.position[2] + intersect_dist * view_vec[2],
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
        let look_vector_raw = subtract(&earth_point, &state.position);
        let look_vector = normalize(&look_vector_raw);

        let (_sat_x, sat_y) =
            self.satellite_body_axes(&state.spin_axis, &state.sun_vec, state.beta_angle);

        // Compute working angles for instrument pointing
        let pixel_angle = {
            let cross_spin_look = cross(&state.spin_axis, &look_vector);
            let cross_y_cross = cross(&sat_y, &cross_spin_look);
            let mut angle = angle_between(&sat_y, &cross_spin_look);
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
