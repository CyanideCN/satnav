use crate::math3d::{add_scaled, change_basis, dot, magnitude, mat3_mul_vec3, Vec3};

const EARTH_ROT_RATE: f64 = 0.26251617;
const EARTH_RAD_PRODUCT: f64 = 40546851.22;
const EARTH_EQ_RADIUS_SQ: f64 = 40683833.48;
const EARTH_POL_RADIUS_SQ: f64 = 40410330.18;
const POL_TO_EQ_RADIUS_RATIO_SQ: f64 = EARTH_POL_RADIUS_SQ / EARTH_EQ_RADIUS_SQ;
const ONE_MINUS_RADIUS_RATIO_SQ: f64 = 1.0 - POL_TO_EQ_RADIUS_RATIO_SQ;
const EARTH_RAD_KM: f64 = 6371.22;
const EARTH_RAD_KM_SQ: f64 = EARTH_RAD_KM * EARTH_RAD_KM;
const PI: f64 = std::f64::consts::PI;
const DEG_TO_RAD: f64 = PI / 180.0;
const STELLAR_HOUR_ANGLE: f64 = 100.26467;
const REFERENCE_DAY: i32 = 74001;
const REFERENCE_TIME: i32 = 0;
const SOLAR_SIDEREAL_RATIO: f64 = 1.00273791;
const GRAVITATIONAL_CONSTANT: f64 = 0.07436574;

fn packed_yymmdd_to_yyddd(packed_date: i32) -> i32 {
    let cumulative_days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];

    let year = (packed_date / 10000) % 100;
    let mut month = (packed_date / 100) % 100;
    let day = packed_date % 100;

    if !(1..=12).contains(&month) {
        month = 1;
    }

    let mut julian_day = day + cumulative_days[(month - 1) as usize];
    if year % 4 == 0 && month > 2 {
        julian_day += 1;
    }

    1000 * year + julian_day
}

fn days_in_year(year: i32) -> i32 {
    366 - (year % 4 + 3) / 4
}

fn decimal_hours_to_packed_hms(decimal_hours: f64) -> i32 {
    let hours = decimal_hours as i32;
    let minutes = ((decimal_hours - hours as f64) * 60.0) as i32;
    let seconds = (((decimal_hours - hours as f64) * 60.0 - minutes as f64) * 60.0) as i32;
    hours * 10000 + minutes * 100 + seconds
}

fn packed_dms_to_decimal(packed_value: i32) -> f64 {
    let absolute_value = packed_value.abs();
    let degrees_or_hours = (absolute_value / 10000) as f64;
    let minutes = ((absolute_value / 100) % 100) as f64 / 60.0;
    let seconds = (absolute_value % 100) as f64 / 3600.0;

    let mut decimal_value = degrees_or_hours + minutes + seconds;
    if packed_value < 0 {
        decimal_value = -decimal_value;
    }

    decimal_value
}

fn yyddd_components(packed_date: i32) -> (i32, i32) {
    ((packed_date / 1000) % 100, packed_date % 1000)
}

fn elapsed_days_since_year_zero(packed_date: i32) -> i32 {
    let (year, day_of_year) = yyddd_components(packed_date);
    let leap_days = (year - 1) / 4 + 1;
    365 * (year - 1) + leap_days + day_of_year - 1
}

fn perigee_epoch(
    epoch_date: i32,
    epoch_time: i32,
    semi_major_axis: f64,
    eccentricity: f64,
    mean_anomaly_deg: f64,
) -> (i32, i32) {
    let earth_radius = 6378.388;
    let mean_motion = GRAVITATIONAL_CONSTANT * (earth_radius / semi_major_axis).sqrt().powi(3);
    let mean_anomaly_rad = DEG_TO_RAD * mean_anomaly_deg;

    let mut adjusted_time = packed_dms_to_decimal(epoch_time)
        - (mean_anomaly_rad - eccentricity * mean_anomaly_rad.sin()) / (60.0 * mean_motion);
    let mut day_offset = 0;

    if adjusted_time > 48.0 {
        adjusted_time -= 48.0;
        day_offset = 2;
    } else if adjusted_time > 24.0 {
        adjusted_time -= 24.0;
        day_offset = 1;
    } else if adjusted_time < -24.0 {
        adjusted_time += 48.0;
        day_offset = -2;
    } else if adjusted_time < 0.0 {
        adjusted_time += 24.0;
        day_offset = -1;
    }

    let updated_epoch_time = decimal_hours_to_packed_hms(adjusted_time);
    if day_offset == 0 {
        return (epoch_date, updated_epoch_time);
    }

    let (mut year, mut day_of_year) = yyddd_components(epoch_date);
    day_of_year += day_offset;

    if day_of_year < 1 {
        year -= 1;
        day_of_year += days_in_year(year);
    } else {
        let days_this_year = days_in_year(year);
        if day_of_year > days_this_year {
            year += 1;
            day_of_year -= days_this_year;
        }
    }

    (1000 * year + day_of_year, updated_epoch_time)
}

fn minutes_between(date1: i32, time1: i32, date2: i32, time2: i32) -> f64 {
    let total_minutes_1 =
        1440.0 * elapsed_days_since_year_zero(date1) as f64 + 60.0 * packed_dms_to_decimal(time1);
    let total_minutes_2 =
        1440.0 * elapsed_days_since_year_zero(date2) as f64 + 60.0 * packed_dms_to_decimal(time2);

    total_minutes_2 - total_minutes_1
}

#[allow(dead_code)]
fn celestial_to_earth_longitude(epoch_date: i32, epoch_time: i32, celestial_longitude: f64) -> f64 {
    let elapsed_minutes = minutes_between(epoch_date, epoch_time, REFERENCE_DAY, REFERENCE_TIME);
    let right_ascension_hour_angle =
        celestial_longitude - STELLAR_HOUR_ANGLE + (elapsed_minutes * SOLAR_SIDEREAL_RATIO) / 4.0;

    let mut earth_longitude = right_ascension_hour_angle % 360.0;
    if earth_longitude < 0.0 {
        earth_longitude += 360.0;
    }

    earth_longitude
}

fn earth_to_celestial_longitude(epoch_date: i32, epoch_time: i32, earth_longitude: f64) -> f64 {
    let elapsed_minutes = minutes_between(REFERENCE_DAY, REFERENCE_TIME, epoch_date, epoch_time);
    let right_ascension_hour_angle =
        earth_longitude + (elapsed_minutes * SOLAR_SIDEREAL_RATIO) / 4.0 + STELLAR_HOUR_ANGLE;

    let mut celestial_longitude = right_ascension_hour_angle % 360.0;
    if celestial_longitude < 0.0 {
        celestial_longitude += 360.0;
    }

    celestial_longitude
}

fn xyz_to_geodetic_latlon(point: &Vec3) -> (f64, f64) {
    let [x, y, z] = *point;
    if x == 0.0 && y == 0.0 && z == 0.0 {
        return (100.0, 200.0);
    }

    let geocentric_angle = (z / (x.powi(2) + y.powi(2)).sqrt()).atan();
    let latitude = ((EARTH_EQ_RADIUS_SQ * geocentric_angle.sin())
        .atan2(EARTH_POL_RADIUS_SQ * geocentric_angle.cos()))
        / DEG_TO_RAD;

    let mut longitude = -(y.atan2(x)) / DEG_TO_RAD;
    if longitude < -180.0 {
        longitude += 360.0;
    } else if longitude > 180.0 {
        longitude -= 360.0;
    }

    (latitude, longitude)
}

fn geodetic_latlon_to_xyz(latitude: f64, longitude: f64) -> Vec3 {
    let geocentric_latitude = (EARTH_POL_RADIUS_SQ * (latitude * DEG_TO_RAD).sin())
        .atan2(EARTH_EQ_RADIUS_SQ * (latitude * DEG_TO_RAD).cos());
    let longitude_rad = -DEG_TO_RAD * longitude;
    let (sin_latitude, cos_latitude) = geocentric_latitude.sin_cos();
    let (sin_longitude, cos_longitude) = longitude_rad.sin_cos();
    let tan_latitude_sq = (sin_latitude / cos_latitude).powi(2);
    let radius = EARTH_RAD_PRODUCT
        * ((1.0 + tan_latitude_sq) / (EARTH_POL_RADIUS_SQ + EARTH_EQ_RADIUS_SQ * tan_latitude_sq))
            .sqrt();

    let x = radius * cos_latitude * cos_longitude;
    let y = radius * cos_latitude * sin_longitude;
    let z = radius * sin_latitude;

    [x, y, z]
}

fn rotate_xy(vector: &Vec3, angle: f64) -> Vec3 {
    let (sin_angle, cos_angle) = angle.sin_cos();
    [
        cos_angle * vector[0] - sin_angle * vector[1],
        sin_angle * vector[0] + cos_angle * vector[1],
        vector[2],
    ]
}

#[derive(Debug)]
pub struct GOESNavigation {
    pub navigation_day: i32,
    pub epoch_date: i32,
    pub epoch_time: i32,
    pub semi_major_axis: f64,
    pub orbital_eccentricity: f64,
    pub orbital_inclination: f64,
    pub mean_anomaly_deg: f64,
    pub argument_of_perigee_deg: f64,
    pub ascending_node_deg: f64,
    pub declination_deg: f64,
    pub right_ascension_deg: f64,
    pub image_center_line: f64,
    pub spin_period: f64,
    pub line_sweep_angle_deg: f64,
    pub encoded_scan_line_count: i32,
    pub element_sweep_angle_deg: f64,
    pub element_count: i32,
    pub pitch_deg: f64,
    pub yaw_deg: f64,
    pub roll_deg: f64,
    pub skew: f64,
    pub element_offset: f64,
    pub element_drift: f64,
    pub scan_start_number: f64,
    pub scan_start_time: f64,
    pub scan_end_number: f64,
    pub scan_end_time: f64,
    pub uses_legacy_timing: bool,
    pub sensor_count: i32,
    pub total_lines: i32,
    pub line_angle_step: f64,
    pub total_elements: i32,
    pub element_angle_step: f64,
    pub image_center_element: f64,
    pub line_rotation_11: f64,
    pub line_rotation_13: f64,
    pub line_rotation_21: f64,
    pub line_rotation_23: f64,
    pub line_rotation_31: f64,
    pub line_rotation_33: f64,
    pub line_rotation_norm_sq: f64,
    pub line_rotation_phase: f64,
    pub frame_time_scale: f64,
    pub celestial_transform_11: f64,
    pub celestial_transform_12: f64,
    pub celestial_transform_13: f64,
    pub celestial_transform_21: f64,
    pub celestial_transform_22: f64,
    pub celestial_transform_23: f64,
    pub celestial_transform_31: f64,
    pub celestial_transform_32: f64,
    pub celestial_transform_33: f64,
    pub earth_reference_rotation: f64,
    pub picture_time: f64,
}

impl GOESNavigation {
    pub fn new(nav_block: &[i32]) -> GOESNavigation {
        let epoch_date = packed_yymmdd_to_yyddd(nav_block[4]);
        let epoch_time =
            100 * (nav_block[5] / 100) + (0.6 * (nav_block[5] % 100) as f64).round() as i32;

        let mut image_center_line = nav_block[14] as f64;
        if image_center_line > 1_000_000.0 {
            image_center_line /= 10_000.0;
        }

        if nav_block[12] == 0 && nav_block[13] == 0 && nav_block[14] == 0 {
            panic!("Invalid ascension/declination parameters");
        }
        if nav_block[15] == 0 {
            panic!("Invalid spin period");
        }

        let mut navigation = GOESNavigation {
            navigation_day: nav_block[1] % 100000,
            epoch_date,
            epoch_time,
            semi_major_axis: nav_block[6] as f64 / 100.0,
            orbital_eccentricity: nav_block[7] as f64 / 1_000_000.0,
            orbital_inclination: nav_block[8] as f64 / 1000.0,
            mean_anomaly_deg: nav_block[9] as f64 / 1000.0,
            argument_of_perigee_deg: nav_block[10] as f64 / 1000.0,
            ascending_node_deg: nav_block[11] as f64 / 1000.0,
            declination_deg: packed_dms_to_decimal(nav_block[12]),
            right_ascension_deg: packed_dms_to_decimal(nav_block[13]),
            image_center_line,
            spin_period: nav_block[15] as f64 / 1000.0,
            line_sweep_angle_deg: packed_dms_to_decimal(nav_block[16]),
            encoded_scan_line_count: nav_block[17],
            element_sweep_angle_deg: packed_dms_to_decimal(nav_block[18]),
            element_count: nav_block[19],
            pitch_deg: packed_dms_to_decimal(nav_block[20]),
            yaw_deg: packed_dms_to_decimal(nav_block[21]),
            roll_deg: packed_dms_to_decimal(nav_block[22]),
            skew: if nav_block[28] == 0x80808080u32 as i32 {
                0.0
            } else {
                nav_block[28] as f64 / 100000.0
            },
            element_offset: nav_block[38] as f64 / 100.0,
            element_drift: nav_block[39] as f64 / 100.0,
            scan_start_number: 0.0,
            scan_start_time: packed_dms_to_decimal(nav_block[31]),
            scan_end_number: 0.0,
            scan_end_time: 0.0,
            uses_legacy_timing: false,
            sensor_count: 1,
            total_lines: 0,
            line_angle_step: 0.0,
            total_elements: 0,
            element_angle_step: 0.0,
            image_center_element: 0.0,
            line_rotation_11: 0.0,
            line_rotation_13: 0.0,
            line_rotation_21: 0.0,
            line_rotation_23: 0.0,
            line_rotation_31: 0.0,
            line_rotation_33: 0.0,
            line_rotation_norm_sq: 0.0,
            line_rotation_phase: 0.0,
            frame_time_scale: 0.0,
            celestial_transform_11: 0.0,
            celestial_transform_12: 0.0,
            celestial_transform_13: 0.0,
            celestial_transform_21: 0.0,
            celestial_transform_22: 0.0,
            celestial_transform_23: 0.0,
            celestial_transform_31: 0.0,
            celestial_transform_32: 0.0,
            celestial_transform_33: 0.0,
            earth_reference_rotation: earth_to_celestial_longitude(nav_block[1], 0, 0.0)
                * DEG_TO_RAD,
            picture_time: packed_dms_to_decimal(nav_block[2]),
        };

        let (perigee_epoch_date, perigee_epoch_time) = perigee_epoch(
            epoch_date,
            epoch_time,
            navigation.semi_major_axis,
            navigation.orbital_eccentricity,
            navigation.mean_anomaly_deg,
        );
        navigation.epoch_date = perigee_epoch_date;
        navigation.epoch_time = perigee_epoch_time;

        if navigation.spin_period < 300.0 {
            navigation.spin_period = 60000.0 / navigation.spin_period;
        }

        navigation.sensor_count = (navigation.encoded_scan_line_count / 100000) % 100;
        if navigation.sensor_count < 1 {
            navigation.sensor_count = 1;
        }

        navigation.total_lines =
            navigation.sensor_count * (navigation.encoded_scan_line_count % 100000);
        navigation.line_angle_step =
            DEG_TO_RAD * navigation.line_sweep_angle_deg / (navigation.total_lines - 1) as f64;
        navigation.total_elements = navigation.element_count;
        navigation.element_angle_step = DEG_TO_RAD * navigation.element_sweep_angle_deg
            / (navigation.total_elements - 1) as f64;
        navigation.image_center_element = (1.0 + navigation.total_elements as f64) / 2.0;

        let pitch_rad = DEG_TO_RAD * navigation.pitch_deg;
        let yaw_rad = DEG_TO_RAD * navigation.yaw_deg;
        let roll_rad = DEG_TO_RAD * navigation.roll_deg;
        let skew_angle = navigation
            .skew
            .atan2(navigation.line_angle_step / navigation.element_angle_step);

        let (sin_pitch, cos_pitch) = pitch_rad.sin_cos();
        let (sin_yaw, cos_yaw) = (yaw_rad - skew_angle).sin_cos();
        let (sin_roll, cos_roll) = roll_rad.sin_cos();

        navigation.line_rotation_11 = cos_roll * cos_pitch;
        navigation.line_rotation_13 = sin_yaw * sin_roll * cos_pitch + cos_yaw * sin_pitch;
        navigation.line_rotation_21 = -sin_roll;
        navigation.line_rotation_23 = sin_yaw * cos_roll;
        navigation.line_rotation_31 = -cos_roll * sin_pitch;
        navigation.line_rotation_33 = cos_yaw * cos_pitch - sin_yaw * sin_roll * sin_pitch;
        navigation.line_rotation_norm_sq =
            navigation.line_rotation_31.powi(2) + navigation.line_rotation_33.powi(2);
        navigation.line_rotation_phase = navigation
            .line_rotation_31
            .atan2(navigation.line_rotation_33);
        navigation.frame_time_scale = navigation.spin_period / 3_600_000.0;

        let declination_rad = navigation.declination_deg * DEG_TO_RAD;
        let (sin_declination, cos_declination) = declination_rad.sin_cos();
        let right_ascension_rad = navigation.right_ascension_deg * DEG_TO_RAD;
        let (sin_right_ascension, cos_right_ascension) = right_ascension_rad.sin_cos();

        navigation.celestial_transform_11 = -sin_right_ascension;
        navigation.celestial_transform_12 = cos_right_ascension;
        navigation.celestial_transform_13 = 0.0;
        navigation.celestial_transform_21 = -sin_declination * cos_right_ascension;
        navigation.celestial_transform_22 = -sin_declination * sin_right_ascension;
        navigation.celestial_transform_23 = cos_declination;
        navigation.celestial_transform_31 = cos_declination * cos_right_ascension;
        navigation.celestial_transform_32 = cos_declination * sin_right_ascension;
        navigation.celestial_transform_33 = sin_declination;

        let satellite_series = nav_block[1] / 100000;
        if (satellite_series > 25 || satellite_series == 12) && nav_block[30] > 0 {
            navigation.scan_start_number = nav_block[30] as f64;
            navigation.scan_start_time = packed_dms_to_decimal(nav_block[31]);
            navigation.scan_end_number = nav_block[34] as f64;
            navigation.scan_end_time = packed_dms_to_decimal(nav_block[35]);
        } else {
            navigation.scan_start_number = 1.0;
            navigation.scan_start_time = 0.0;
            navigation.scan_end_number = 2.0;
            navigation.scan_end_time = 0.0;
            navigation.uses_legacy_timing = true;
        }

        navigation
    }

    fn satellite_vector(&self, satellite_time: f64) -> Vec3 {
        let earth_radius = 6378.388;

        let orbital_inclination_rad = DEG_TO_RAD * self.orbital_inclination;
        let argument_of_perigee_rad = DEG_TO_RAD * self.argument_of_perigee_deg;
        let ascending_node_rad = DEG_TO_RAD * self.ascending_node_deg;

        let (sin_inclination, cos_inclination) = orbital_inclination_rad.sin_cos();
        let (mut sin_perigee, mut cos_perigee) = argument_of_perigee_rad.sin_cos();
        sin_perigee *= self.semi_major_axis;
        cos_perigee *= self.semi_major_axis;
        let (sin_ascending_node, cos_ascending_node) = ascending_node_rad.sin_cos();

        let position_basis_x =
            cos_perigee * cos_ascending_node - sin_perigee * sin_ascending_node * cos_inclination;
        let position_basis_y =
            cos_perigee * sin_ascending_node + sin_perigee * cos_ascending_node * cos_inclination;
        let position_basis_z = sin_perigee * sin_inclination;

        let velocity_basis_x =
            -sin_perigee * cos_ascending_node - cos_perigee * sin_ascending_node * cos_inclination;
        let velocity_basis_y =
            -sin_perigee * sin_ascending_node + cos_perigee * cos_ascending_node * cos_inclination;
        let velocity_basis_z = cos_perigee * sin_inclination;

        let orbital_eccentricity_scale =
            ((1.0 - self.orbital_eccentricity) * (1.0 + self.orbital_eccentricity)).sqrt();
        let mean_motion =
            GRAVITATIONAL_CONSTANT * earth_radius * (earth_radius / self.semi_major_axis).sqrt()
                / self.semi_major_axis;

        let epoch_elapsed_minutes = 1440.0 * elapsed_days_since_year_zero(self.epoch_date) as f64
            + 60.0 * packed_dms_to_decimal(self.epoch_time);
        let navigation_elapsed_minutes =
            elapsed_days_since_year_zero(self.navigation_day) as f64 * 1440.0;
        let elapsed_minutes =
            navigation_elapsed_minutes - epoch_elapsed_minutes + satellite_time * 60.0;

        let mean_anomaly = mean_motion * elapsed_minutes;
        let mut previous_eccentric_anomaly = mean_anomaly;
        let mut eccentric_anomaly =
            mean_anomaly + self.orbital_eccentricity * previous_eccentric_anomaly.sin();

        for _ in 0..20 {
            eccentric_anomaly =
                mean_anomaly + self.orbital_eccentricity * previous_eccentric_anomaly.sin();
            if (eccentric_anomaly - previous_eccentric_anomaly).abs() < 1.0e-8 {
                break;
            }
            previous_eccentric_anomaly = eccentric_anomaly;
        }

        let (sin_eccentric_anomaly, cos_eccentric_anomaly) = eccentric_anomaly.sin_cos();
        let orbit_x = cos_eccentric_anomaly - self.orbital_eccentricity;
        let orbit_y = orbital_eccentricity_scale * sin_eccentric_anomaly;

        let z = orbit_x * position_basis_z + orbit_y * velocity_basis_z;
        let y = orbit_x * position_basis_y + orbit_y * velocity_basis_y;
        let x = orbit_x * position_basis_x + orbit_y * velocity_basis_x;

        [x, y, z]
    }

    pub fn image_to_geographic(&self, line: f64, pixel: f64) -> (f64, f64) {
        let nearest_line = (line + 0.5).floor();
        let frame_index = (nearest_line - 1.0) / self.sensor_count as f64 + 1.0;
        let frame_time = self.frame_time_scale * frame_index;
        let sample_time = frame_time + self.picture_time;

        let satellite_position = self.satellite_vector(sample_time);

        let line_angle = (line - self.image_center_line) * self.line_angle_step;
        let element_angle = (pixel - self.image_center_element
            + self.element_offset
            + self.element_drift * sample_time)
            * self.element_angle_step;

        let celestial_transform = [
            [
                self.celestial_transform_11,
                self.celestial_transform_12,
                self.celestial_transform_13,
            ],
            [
                self.celestial_transform_21,
                self.celestial_transform_22,
                self.celestial_transform_23,
            ],
            [
                self.celestial_transform_31,
                self.celestial_transform_32,
                self.celestial_transform_33,
            ],
        ];
        let satellite_celestial = mat3_mul_vec3(&celestial_transform, &satellite_position);
        let earth_rotation_offset = satellite_celestial[1].atan2(satellite_celestial[0]) + PI;
        let adjusted_element_angle = element_angle - earth_rotation_offset;

        let (sin_line, cos_line) = line_angle.sin_cos();
        let forward_line_rotation = [
            [self.line_rotation_11, 0.0, -self.line_rotation_13],
            [self.line_rotation_21, 0.0, -self.line_rotation_23],
            [self.line_rotation_31, 0.0, -self.line_rotation_33],
        ];
        let line_direction = mat3_mul_vec3(&forward_line_rotation, &[cos_line, 0.0, sin_line]);
        let scan_direction = rotate_xy(&line_direction, -adjusted_element_angle);
        let earth_ray = change_basis(
            &celestial_transform[0],
            &celestial_transform[1],
            &celestial_transform[2],
            &scan_direction,
        );

        let quadratic_a =
            POL_TO_EQ_RADIUS_RATIO_SQ + ONE_MINUS_RADIUS_RATIO_SQ * earth_ray[2] * earth_ray[2];
        let quadratic_b = 2.0
            * ((earth_ray[0] * satellite_position[0] + earth_ray[1] * satellite_position[1])
                * POL_TO_EQ_RADIUS_RATIO_SQ
                + earth_ray[2] * satellite_position[2]);
        let quadratic_c = (satellite_position[0] * satellite_position[0]
            + satellite_position[1] * satellite_position[1])
            * POL_TO_EQ_RADIUS_RATIO_SQ
            + satellite_position[2] * satellite_position[2]
            - EARTH_POL_RADIUS_SQ;
        let discriminant = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c;

        if discriminant < 1.0 {
            return (f64::NAN, f64::NAN);
        }

        let distance_along_ray = -(quadratic_b + discriminant.sqrt()) / (2.0 * quadratic_a);
        let earth_intersection =
            add_scaled(&satellite_position, 1.0, &earth_ray, distance_along_ray);

        let rotated_earth = rotate_xy(
            &earth_intersection,
            -(EARTH_ROT_RATE * sample_time + self.earth_reference_rotation),
        );

        let (latitude, longitude_west_positive) = xyz_to_geodetic_latlon(&rotated_earth);
        (latitude, -longitude_west_positive)
    }

    pub fn get_subpoint(&self) -> (f64, f64) {
        let sample_time = self.scan_start_time;
        let satellite_position = self.satellite_vector(sample_time);
        let rotated_position = rotate_xy(
            &satellite_position,
            -(EARTH_ROT_RATE * sample_time + self.earth_reference_rotation),
        );
        let (latitude, longitude_west_positive) = xyz_to_geodetic_latlon(&rotated_position);
        (latitude, -longitude_west_positive)
    }

    pub fn geographic_to_image(&self, lat: f64, lon: f64) -> (f64, f64) {
        let west_positive_lon = -lon;
        let earth_position = geodetic_latlon_to_xyz(lat, west_positive_lon);
        let celestial_transform = [
            [
                self.celestial_transform_11,
                self.celestial_transform_12,
                self.celestial_transform_13,
            ],
            [
                self.celestial_transform_21,
                self.celestial_transform_22,
                self.celestial_transform_23,
            ],
            [
                self.celestial_transform_31,
                self.celestial_transform_32,
                self.celestial_transform_33,
            ],
        ];

        let mut estimated_line = f64::NAN;
        let mut first_pass_line = 910.0;
        let mut cached_satellite_time = -99999.0;
        let mut satellite_position = [0.0, 0.0, 0.0];
        let mut rotated_earth = [0.0, 0.0, 0.0];
        let mut satellite_distance = 0.0;
        let mut line_of_sight_length = 0.0;
        let mut sample_time = self.scan_start_time;

        for pass in 0..2 {
            if (sample_time - cached_satellite_time).abs() >= 0.0005 {
                satellite_position = self.satellite_vector(sample_time);
                cached_satellite_time = sample_time;
                satellite_distance = magnitude(&satellite_position);
            }

            rotated_earth = rotate_xy(
                &earth_position,
                EARTH_ROT_RATE * sample_time + self.earth_reference_rotation,
            );

            let line_of_sight = [
                rotated_earth[0] - satellite_position[0],
                rotated_earth[1] - satellite_position[1],
                rotated_earth[2] - satellite_position[2],
            ];
            let transformed_line_of_sight = mat3_mul_vec3(&celestial_transform, &line_of_sight);
            line_of_sight_length = magnitude(&line_of_sight);

            let normalized_line_of_sight_z = transformed_line_of_sight[2] / line_of_sight_length;
            let scan_angle = normalized_line_of_sight_z
                .atan2((self.line_rotation_norm_sq - normalized_line_of_sight_z.powi(2)).sqrt())
                - self.line_rotation_phase;
            estimated_line = self.image_center_line - scan_angle / self.line_angle_step;

            if pass == 0 {
                sample_time = self.scan_end_time;
                first_pass_line = estimated_line;
            }
        }

        let scan_number =
            ((first_pass_line + estimated_line) / 2.0 - 1.0) / self.sensor_count as f64;
        let scan_fraction = (scan_number - self.scan_start_number)
            / (self.scan_end_number - self.scan_start_number);
        estimated_line = first_pass_line + scan_fraction * (estimated_line - first_pass_line);
        sample_time =
            self.scan_start_time + self.frame_time_scale * (scan_number - self.scan_start_number);

        satellite_position = self.satellite_vector(sample_time);

        let visibility_dot = dot(&rotated_earth, &satellite_position);
        let horizon_threshold = 0.0001 * EARTH_RAD_KM * satellite_distance + EARTH_RAD_KM_SQ;
        let mut estimated_element = f64::NAN;

        if visibility_dot >= horizon_threshold {
            let satellite_celestial = mat3_mul_vec3(&celestial_transform, &satellite_position);
            rotated_earth = rotate_xy(
                &earth_position,
                EARTH_ROT_RATE * sample_time + self.earth_reference_rotation,
            );

            let line_of_sight = [
                rotated_earth[0] - satellite_position[0],
                rotated_earth[1] - satellite_position[1],
                rotated_earth[2] - satellite_position[2],
            ];
            let transformed_line_of_sight = mat3_mul_vec3(&celestial_transform, &line_of_sight);

            let perpendicular_line_of_sight =
                (line_of_sight_length.powi(2) - transformed_line_of_sight[2].powi(2)).sqrt();
            let satellite_xy_length =
                magnitude(&[satellite_celestial[0], satellite_celestial[1], 0.0]);
            line_of_sight_length = magnitude(&line_of_sight);

            let normalized_line_of_sight_z = transformed_line_of_sight[2] / line_of_sight_length;
            let line_angle = normalized_line_of_sight_z
                .atan2((self.line_rotation_norm_sq - normalized_line_of_sight_z.powi(2)).sqrt())
                - self.line_rotation_phase;
            let (sin_line, cos_line) = line_angle.sin_cos();

            let inverse_line_rotation = [
                [self.line_rotation_11, 0.0, self.line_rotation_13],
                [self.line_rotation_21, 0.0, self.line_rotation_23],
                [self.line_rotation_31, 0.0, self.line_rotation_33],
            ];
            let inverse_line_direction =
                mat3_mul_vec3(&inverse_line_rotation, &[cos_line, 0.0, sin_line]);

            estimated_element = self.image_center_element
                + ((satellite_celestial[0] * transformed_line_of_sight[1]
                    - satellite_celestial[1] * transformed_line_of_sight[0])
                    / (perpendicular_line_of_sight * satellite_xy_length))
                    .asin()
                    / self.element_angle_step;
            estimated_element += inverse_line_direction[1].atan2(inverse_line_direction[0])
                / self.element_angle_step;
            estimated_element -= self.element_offset + self.element_drift * sample_time;
        }

        (estimated_line, estimated_element)
    }
}
