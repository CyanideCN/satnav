const EMEGA: f64 = 0.26251617;
const AB: f64 = 40546851.22;
const ASQ: f64 = 40683833.48;
const BSQ: f64 = 40410330.18;
const BASQ: f64 = BSQ / ASQ;
const ONEMSQ: f64 = 1.0 - BASQ;
const RE: f64 = 6371.22;
const RSQ: f64 = RE * RE;
const PI: f64 = std::f64::consts::PI;
const RDPDG: f64 = PI / 180.;
const SHA: f64 = 100.26467; // Stellar Hour Angle
const IRAYD: i32 = 74001; // Reference year and day (YYDDD format)
const IRAHMS: i32 = 0; // Reference time (HHMMSS format)
const SOLSID: f64 = 1.00273791; // Solar sidereal time ratio
const GRAVCON: f64 = 0.07436574;

fn icon1(date: i32) -> i32 {
    // Define the number of days cumulatively at the start of each month
    let cumulative_days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];

    // Extract year, month, and day from the input date (in YYMMDD format)
    let year = (date / 10000) % 100;
    let mut month = (date / 100) % 100;
    let day = date % 100;

    // Validate month (ensure month is between 1 and 12)
    if month < 1 || month > 12 {
        month = 1;
    }

    // Calculate the Julian day for the given month and day
    let mut julian_day = day + cumulative_days[(month - 1) as usize];

    // Check if the year is a leap year and adjust the Julian day if needed
    if year % 4 == 0 && month > 2 {
        julian_day += 1;
    }

    // Combine year and Julian day
    1000 * year + julian_day
}

fn is_leap_year(year: i32) -> i32 {
    // Determine if the given year is a leap year (returns 366 if leap year, 365 otherwise)
    366 - (year % 4 + 3) / 4
}

fn float_time_to_int(time: f64) -> i32 {
    // Convert time in HH.MMSS format to an integer in the format HHMMSS
    let hours = time as i32;
    let minutes = ((time - hours as f64) * 60.0) as i32;
    let seconds = (((time - hours as f64) * 60.0 - minutes as f64) * 60.0) as i32;
    hours * 10000 + minutes * 100 + seconds
}

fn flalo(value: i32) -> f64 {
    // Convert an integer angle (DDDMMSS) or time (HHMMSS) to a floating point value
    let abs_value = value.abs();

    // Convert the integer to degrees/hours, minutes, and seconds
    let degrees_or_hours = (abs_value / 10000) as f64;
    let minutes = ((abs_value / 100) % 100) as f64 / 60.0;
    let seconds = (abs_value % 100) as f64 / 3600.0;

    // Combine the parts into a floating point result
    let mut result = degrees_or_hours + minutes + seconds;

    // Apply the sign of the original value
    if value < 0 {
        result = -result;
    }

    result
}

fn perigee_time(
    epoch_year: i32,
    epoch_time: i32,
    semi_major_axis: f64,
    eccentricity: f64,
    mean_anomaly_deg: f64,
) -> (i32, i32) {
    // Constants
    let earth_radius = 6378.388;
    // Calculate mean motion (xmmc)
    let mean_motion = GRAVCON * (earth_radius / semi_major_axis).sqrt().powi(3);
    // Convert mean anomaly from degrees to radians
    let mean_anomaly_rad = RDPDG * mean_anomaly_deg;
    // Calculate time in minutes (Kepler's equation approximation)
    let mut time_offset =
        (mean_anomaly_rad - eccentricity * mean_anomaly_rad.sin()) / (60.0 * mean_motion);
    // Convert input time (HHMMSS format) to float and subtract time offset
    let epoch_time_float = flalo(epoch_time);
    time_offset = epoch_time_float - time_offset;
    let mut day_offset = 0;
    // Adjust the time based on the day offset
    if time_offset > 48.0 {
        time_offset -= 48.0;
        day_offset = 2;
    } else if time_offset > 24.0 {
        time_offset -= 24.0;
        day_offset = 1;
    } else if time_offset < -24.0 {
        time_offset += 48.0;
        day_offset = -2;
    } else if time_offset < 0.0 {
        time_offset += 24.0;
        day_offset = -1;
    }
    // Convert the adjusted time back to integer format (HHMMSS)
    let updated_epoch_time = float_time_to_int(time_offset);
    // Handle day adjustments
    if day_offset == 0 {
        return (epoch_year, updated_epoch_time);
    }
    // Extract year and day from epoch year (YYDDD format)
    let mut year = (epoch_year / 1000) % 100;
    let mut day_of_year = epoch_year % 1000;
    // Adjust day of the year based on the day offset
    day_of_year += day_offset;
    if day_of_year < 1 {
        year -= 1;
        day_of_year += is_leap_year(year);
    } else {
        let total_days_in_year = is_leap_year(year);
        if day_of_year > total_days_in_year {
            year += 1;
            day_of_year -= total_days_in_year;
        }
    }
    // Recombine the year and day into the new epoch year (YYDDD format)
    let updated_epoch_year = 1000 * year + day_of_year;

    (updated_epoch_year, updated_epoch_time)
}

fn time_difference(date1: i32, time1: i32, date2: i32, time2: i32) -> f64 {
    // Extract year and day for the first date (YYDDD format)
    let year1 = (date1 / 1000) % 100;
    let day1 = date1 % 1000;
    let leap_factor1 = (year1 - 1) / 4 + 1;
    // Calculate total number of days since year 0 for the first date
    let total_days1 = 365 * (year1 - 1) + leap_factor1 + day1 - 1;
    // Extract year and day for the second date (YYDDD format)
    let year2 = (date2 / 1000) % 100;
    let day2 = date2 % 1000;
    let leap_factor2 = (year2 - 1) / 4 + 1;
    // Calculate total number of days since year 0 for the second date
    let total_days2 = 365 * (year2 - 1) + leap_factor2 + day2 - 1;
    // Convert times (HHMMSS format) to float values (in minutes)
    let time_in_minutes1 = 1440.0 * total_days1 as f64 + 60.0 * flalo(time1);
    let time_in_minutes2 = 1440.0 * total_days2 as f64 + 60.0 * flalo(time2);
    // Calculate and return the difference in minutes

    time_in_minutes2 - time_in_minutes1
}

#[allow(dead_code)]
fn earth_longitude(epoch_date: i32, epoch_time: i32, celestial_longitude: f64) -> f64 {
    // Calculate the Right Ascension Hour Angle (RAHA)
    let time_difference_minutes = time_difference(epoch_date, epoch_time, IRAYD, IRAHMS);
    let right_ascension_hour_angle =
        celestial_longitude - SHA + (time_difference_minutes * SOLSID) / 4.0;
    // Ensure the Earth Longitude (RAE) is within 0 to 360 degrees
    let mut earth_longitude = right_ascension_hour_angle % 360.0;
    if earth_longitude < 0.0 {
        earth_longitude += 360.0;
    }

    earth_longitude
}

fn celestial_longitude(epoch_date: i32, epoch_time: i32, earth_longitude: f64) -> f64 {
    // Calculate the Right Ascension Hour Angle (RAHA)
    let time_difference_minutes = time_difference(IRAYD, IRAHMS, epoch_date, epoch_time);
    let right_ascension_hour_angle =
        earth_longitude + (time_difference_minutes * SOLSID) / 4.0 + SHA;
    // Ensure the Celestial Longitude (RAC) is within 0 to 360 degrees
    let mut celestial_longitude = right_ascension_hour_angle % 360.0;
    if celestial_longitude < 0.0 {
        celestial_longitude += 360.0;
    }

    celestial_longitude
}

fn xyz_to_latlon(x: f64, y: f64, z: f64) -> (f64, f64) {
    // Handle the case where all coordinates are zero
    if x == 0.0 && y == 0.0 && z == 0.0 {
        return (100.0, 200.0); // Default values (as in the original Fortran code)
    }
    // Calculate the geodetic latitude
    let angle = (z / (x.powi(2) + y.powi(2)).sqrt()).atan();
    let latitude = ((ASQ * angle.sin()).atan2(BSQ * angle.cos())) / RDPDG;
    // Calculate the longitude
    let mut longitude = -(y.atan2(x)) / RDPDG;
    // Adjust longitude to be in the range of -180 to 180 degrees
    if longitude < -180.0 {
        longitude += 360.0;
    } else if longitude > 180.0 {
        longitude -= 360.0;
    }

    (latitude, longitude)
}

fn latlon_to_xyz(xlat: f64, xlon: f64) -> (f64, f64, f64) {
    let ylat = (BSQ * (xlat * RDPDG).sin()).atan2(ASQ * (xlat * RDPDG).cos());
    let ylon = -RDPDG * xlon;
    let (snlt, cslt) = ylat.sin_cos();
    let (snln, csln) = ylon.sin_cos();
    let tnlt = (snlt / cslt).powi(2);
    let r = AB * ((1.0 + tnlt) / (BSQ + ASQ * tnlt)).sqrt();
    let x = r * cslt * csln;
    let y = r * cslt * snln;
    let z = r * snlt;

    (x, y, z)
}

#[derive(Debug)]
pub struct GOESNavigation {
    pub nav_day: i32,
    pub epoch_date: i32,
    pub epoch_hour: i32,
    pub semi_major_axis: f64,
    pub orbital_eccentricity: f64,
    pub orbital_inclination: f64,
    pub mean_anomaly: f64,
    pub perigee: f64,
    pub ascending_node: f64,
    pub declination: f64,
    pub right_ascension: f64,
    pub pic_center_line_num: f64,
    pub spin_period: f64,
    pub line_swp_angle: f64,
    pub scan_line_num: i32,
    pub ele_swp_angle: f64,
    pub ele_num: i32,
    pub pitch: f64,
    pub yaw: f64,
    pub roll: f64,
    pub skew: f64,
    pub gamma: f64,
    pub gamma_dot: f64,
    pub scan1: f64,
    pub time1: f64,
    pub scan2: f64,
    pub time2: f64,
    pub iold: i32,
    pub numsen: i32,
    pub total_lines: i32,
    pub radlin: f64,
    pub total_element: i32,
    pub radele: f64,
    pub picele: f64,
    pub rotm11: f64,
    pub rotm13: f64,
    pub rotm21: f64,
    pub rotm23: f64,
    pub rotm31: f64,
    pub rotm33: f64,
    pub rfact: f64,
    pub roasin: f64,
    pub tmpscl: f64,
    pub b11: f64,
    pub b12: f64,
    pub b13: f64,
    pub b21: f64,
    pub b22: f64,
    pub b23: f64,
    pub b31: f64,
    pub b32: f64,
    pub b33: f64,
    pub xref: f64,
    pub picture_time: f64,
}

impl GOESNavigation {
    pub fn new(navblock: &[i32]) -> GOESNavigation {
        // Local variables used for initialization
        let epoch_date = icon1(navblock[4]);
        let epoch_hour =
            100 * (navblock[5] / 100) + (0.6 * (navblock[5] % 100) as f64).round() as i32;
        let mut pic_center_line_num = navblock[14] as f64;
        if pic_center_line_num > 1000000.0 {
            pic_center_line_num /= 10000.0;
        }
        if navblock[12] == 0 && navblock[13] == 0 && navblock[14] == 0 {
            panic!("Invalid ascension/declination parameters");
        }
        if navblock[15] == 0 {
            panic!("Invalid spin period");
        }

        // Initialize the Navigation struct fields
        let mut navigation = GOESNavigation {
            // INTIALIZE NAVCOM
            nav_day: navblock[1] % 100000,
            epoch_date,
            epoch_hour,
            semi_major_axis: navblock[6] as f64 / 100.0,
            orbital_eccentricity: navblock[7] as f64 / 1000000.0,
            orbital_inclination: navblock[8] as f64 / 1000.0,
            mean_anomaly: navblock[9] as f64 / 1000.0,
            perigee: navblock[10] as f64 / 1000.0,
            ascending_node: navblock[11] as f64 / 1000.0,
            declination: flalo(navblock[12]),
            right_ascension: flalo(navblock[13]),
            pic_center_line_num,
            spin_period: navblock[15] as f64 / 1000.0,
            line_swp_angle: flalo(navblock[16]),
            scan_line_num: navblock[17],
            ele_swp_angle: flalo(navblock[18]),
            ele_num: navblock[19],
            pitch: flalo(navblock[20]),
            yaw: flalo(navblock[21]),
            roll: flalo(navblock[22]),
            skew: if navblock[28] == 0x80808080u32 as i32 {
                0.0
            } else {
                navblock[28] as f64 / 100000.0
            },
            gamma: navblock[38] as f64 / 100.0,
            gamma_dot: navblock[39] as f64 / 100.0,
            scan1: 0.0, // Initialized later
            time1: flalo(navblock[31]),
            scan2: 0.0, // Initialized later
            time2: 0.0, // Initialized later
            iold: 0,    // Initialized later
            numsen: 1,  // This will be updated below
            // Note: `numsen` is set to int32 here, which leads to slightly different
            // behavior compared to the python version.
            total_lines: 0,
            radlin: 0.0,
            total_element: 0,
            radele: 0.0,
            picele: 0.0,
            rotm11: 0.0,
            rotm13: 0.0,
            rotm21: 0.0,
            rotm23: 0.0,
            rotm31: 0.0,
            rotm33: 0.0,
            rfact: 0.0,
            roasin: 0.0,
            tmpscl: 0.0,
            b11: 0.0,
            b12: 0.0,
            b13: 0.0,
            b21: 0.0,
            b22: 0.0,
            b23: 0.0,
            b31: 0.0,
            b32: 0.0,
            b33: 0.0,
            xref: celestial_longitude(navblock[1], 0, 0.0) * RDPDG,
            picture_time: flalo(navblock[2]),
        };

        // EPOCH adjustment
        let (new_epoch_date, new_epoch_hour) = perigee_time(
            epoch_date,
            epoch_hour,
            navigation.semi_major_axis,
            navigation.orbital_eccentricity,
            navigation.mean_anomaly,
        );
        navigation.epoch_date = new_epoch_date;
        navigation.epoch_hour = new_epoch_hour;

        if navigation.spin_period < 300.0 {
            navigation.spin_period = 60000.0 / navigation.spin_period;
        }

        navigation.numsen = (navigation.scan_line_num / 100000) % 100;
        if navigation.numsen < 1 {
            navigation.numsen = 1;
        }

        navigation.total_lines = navigation.numsen * (navigation.scan_line_num % 100000);
        navigation.radlin = RDPDG * navigation.line_swp_angle / (navigation.total_lines - 1) as f64;
        navigation.total_element = navigation.ele_num;
        navigation.radele =
            RDPDG * navigation.ele_swp_angle / (navigation.total_element - 1) as f64;
        navigation.picele = (1.0 + navigation.total_element as f64) / 2.0;

        // Compute rotation matrices
        let cpitch = RDPDG * navigation.pitch;
        let cyaw = RDPDG * navigation.yaw;
        let croll = RDPDG * navigation.roll;
        let pskew = navigation.skew.atan2(navigation.radlin / navigation.radele);

        let (sin_pitch, cos_pitch) = cpitch.sin_cos();
        let (sin_yaw, cos_yaw) = (cyaw - pskew).sin_cos();
        let (sin_roll, cos_roll) = croll.sin_cos();

        navigation.rotm11 = cos_roll * cos_pitch;
        navigation.rotm13 = sin_yaw * sin_roll * cos_pitch + cos_yaw * sin_pitch;
        navigation.rotm21 = -sin_roll;
        navigation.rotm23 = sin_yaw * cos_roll;
        navigation.rotm31 = -cos_roll * sin_pitch;
        navigation.rotm33 = cos_yaw * cos_pitch - sin_yaw * sin_roll * sin_pitch;
        navigation.rfact = navigation.rotm31.powi(2) + navigation.rotm33.powi(2);
        navigation.roasin = navigation.rotm31.atan2(navigation.rotm33);
        navigation.tmpscl = navigation.spin_period / 3600000.0;

        // Set celestial sphere transformation matrix based on right ascension and declination
        let dec = navigation.declination * RDPDG;
        let (sindec, cosdec) = dec.sin_cos();
        let ras = navigation.right_ascension * RDPDG;
        let (sinras, cosras) = ras.sin_cos();

        navigation.b11 = -sinras;
        navigation.b12 = cosras;
        navigation.b13 = 0.0;
        navigation.b21 = -sindec * cosras;
        navigation.b22 = -sindec * sinras;
        navigation.b23 = cosdec;
        navigation.b31 = cosdec * cosras;
        navigation.b32 = cosdec * sinras;
        navigation.b33 = sindec;

        // TIME SPECIFIC handling
        let iss = navblock[1] / 100000;
        if (iss > 25 || iss == 12) && navblock[30] > 0 {
            // VAS BIRDS and GMS section
            navigation.scan1 = navblock[30] as f64;
            navigation.time1 = flalo(navblock[31]);
            navigation.scan2 = navblock[34] as f64;
            navigation.time2 = flalo(navblock[35]);
        } else {
            // OLD GOES BIRDS section
            navigation.scan1 = 1.0;
            navigation.time1 = 0.0;
            navigation.scan2 = 2.0;
            navigation.time2 = 0.0;
            navigation.iold = 1;
        }

        navigation
    }

    fn satvec(&self, satellite_time: f64) -> (f64, f64, f64) {
        // Constants
        let earth_radius = 6378.388;

        // Orbital parameters
        let orb_incl_rad = RDPDG * self.orbital_inclination;
        let peri_rad = RDPDG * self.perigee;
        let asc_node_rad = RDPDG * self.ascending_node;

        let (sin_incl, cos_incl) = orb_incl_rad.sin_cos();
        let (mut sin_peri, mut cos_peri) = peri_rad.sin_cos();
        sin_peri *= self.semi_major_axis;
        cos_peri *= self.semi_major_axis;
        let (sin_asnode, cos_asnode) = asc_node_rad.sin_cos();

        let position_x = cos_peri * cos_asnode - sin_peri * sin_asnode * cos_incl;
        let position_y = cos_peri * sin_asnode + sin_peri * cos_asnode * cos_incl;
        let position_z = sin_peri * sin_incl;

        let velocity_x = -sin_peri * cos_asnode - cos_peri * sin_asnode * cos_incl;
        let velocity_y = -sin_peri * sin_asnode + cos_peri * cos_asnode * cos_incl;
        let velocity_z = cos_peri * sin_incl;

        let orbital_eccentricity_sqrt =
            ((1.0 - self.orbital_eccentricity) * (1.0 + self.orbital_eccentricity)).sqrt();
        let mean_motion = GRAVCON * earth_radius * (earth_radius / self.semi_major_axis).sqrt()
            / self.semi_major_axis;

        // Date calculations
        let year = (self.epoch_date / 1000) % 100;
        let day_of_year = self.epoch_date % 1000;
        let leap_year_adj = (year - 1) / 4 + 1;
        let elapsed_days = 365 * (year - 1) + leap_year_adj + day_of_year - 1;
        let elapsed_minutes = 1440.0 * elapsed_days as f64 + 60.0 * flalo(self.epoch_hour);

        let year_nav = (self.nav_day / 1000) % 100;
        let day_nav = self.nav_day % 1000;
        let leap_year_adj_nav = (year_nav - 1) / 4 + 1;
        let elapsed_days_nav = 365 * (year_nav - 1) + leap_year_adj_nav + day_nav - 1;
        let time_diff = elapsed_days_nav as f64 * 1440.0 - elapsed_minutes;

        let sampled_time = satellite_time * 60.0;
        let total_time_diff = time_diff + sampled_time;

        let mean_anomaly = mean_motion * total_time_diff;
        let mut ecc_anom_prev = mean_anomaly;

        // Iterative solution for the Eccentric Anomaly
        let eps = 1.0e-8;
        let mut ecc_anom = mean_anomaly + self.orbital_eccentricity * ecc_anom_prev.sin();
        for _ in 0..20 {
            ecc_anom = mean_anomaly + self.orbital_eccentricity * ecc_anom_prev.sin();
            if (ecc_anom - ecc_anom_prev).abs() < eps {
                break;
            }
            ecc_anom_prev = ecc_anom;
        }

        let (ecc_anom_sin, ecc_anom_cos) = ecc_anom.sin_cos();
        let x_omega = ecc_anom_cos - self.orbital_eccentricity;
        let y_omega = orbital_eccentricity_sqrt * ecc_anom_sin;

        let z = x_omega * position_z + y_omega * velocity_z;
        let y = x_omega * position_y + y_omega * velocity_y;
        let x = x_omega * position_x + y_omega * velocity_x;

        (x, y, z)
    }

    pub fn image_to_geographic(&self, line: f64, pixel: f64) -> (f64, f64) {
        // Convert xlin to nearest integer and calculate related parameters
        let ilin = (line + 0.5).floor();
        let parlin = (ilin - 1.) / self.numsen as f64 + 1.0;
        let framet = self.tmpscl * parlin;
        let samtim = framet + self.picture_time;

        // Call SATVEC to get satellite coordinates
        let (xsat, ysat, zsat) = self.satvec(samtim);

        // Calculations for line and element
        let ylin = (line - self.pic_center_line_num) * self.radlin;
        let yele = (pixel - self.picele + self.gamma + self.gamma_dot * samtim) * self.radele;

        // Coordinate transformations
        let xcor = self.b11 * xsat + self.b12 * ysat + self.b13 * zsat;
        let ycor = self.b21 * xsat + self.b22 * ysat + self.b23 * zsat;
        let rot = (ycor.atan2(xcor)) + PI;
        let yele = yele - rot;

        let (sinlin, coslin) = ylin.sin_cos();
        let (sinele, cosele) = yele.sin_cos();

        let eli = self.rotm11 * coslin - self.rotm13 * sinlin;
        let emi = self.rotm21 * coslin - self.rotm23 * sinlin;
        let eni = self.rotm31 * coslin - self.rotm33 * sinlin;

        let temp = eli;
        let eli = cosele * eli + sinele * emi;
        let emi = -sinele * temp + cosele * emi;

        let elo = self.b11 * eli + self.b21 * emi + self.b31 * eni;
        let emo = self.b12 * eli + self.b22 * emi + self.b32 * eni;
        let eno = self.b13 * eli + self.b23 * emi + self.b33 * eni;

        let aq = BASQ + ONEMSQ * eno * eno;
        let bq = 2.0 * ((elo * xsat + emo * ysat) * BASQ + eno * zsat);
        let cq = (xsat * xsat + ysat * ysat) * BASQ + zsat * zsat - BSQ;
        let rad = bq * bq - 4.0 * aq * cq;

        if rad < 1.0 {
            return (f64::NAN, f64::NAN); // Indicating an error or "off of Earth"
        }

        let s = -(bq + rad.sqrt()) / (2.0 * aq);
        let x = xsat + elo * s;
        let y = ysat + emo * s;
        let z = zsat + eno * s;

        let (st, ct) = (EMEGA * samtim + self.xref).sin_cos();

        let x1 = ct * x + st * y;
        let y1 = -st * x + ct * y;

        let (xpar, ypar) = xyz_to_latlon(x1, y1, z);
        // Return values
        (xpar, -ypar) // Negate ypar if isEastPositive is true
    }

    pub fn get_subpoint(&self) -> (f64, f64) {
        let samtim = self.time1;
        let (x, y, z) = self.satvec(samtim);
        let angle = EMEGA * samtim + self.xref;
        let (st, ct) = angle.sin_cos();
        let x1 = ct * x + st * y;
        let y1 = -st * x + ct * y;
        let (llat, llon) = xyz_to_latlon(x1, y1, z);
        (llat, -llon)
    }

    pub fn geographic_to_image(&self, lat: f64, lon: f64) -> (f64, f64) {
        let mut lon = lon;
        lon *= -1.0; // Assuming east-positive logic as in the original code

        let mut xlin = f64::NAN;
        let mut oldlin = 910.0;
        let mut orbtim = -99999.0;
        let mut xsat = 0.0;
        let mut ysat = 0.0;
        let mut zsat = 0.0;
        let mut x = 0.0;
        let mut y = 0.0;
        let mut xht = 0.0;
        let mut znorm = 0.0;
        let (x1, y1, z) = latlon_to_xyz(lat, lon);

        let mut samtim = self.time1;

        for i in 0..2 {
            if (samtim - orbtim).abs() >= 0.0005 {
                let (new_xsat, new_ysat, new_zsat) = self.satvec(samtim);
                xsat = new_xsat;
                ysat = new_ysat;
                zsat = new_zsat;
                orbtim = samtim;
                xht = (xsat.powi(2) + ysat.powi(2) + zsat.powi(2)).sqrt();
            }

            let (st, ct) = (EMEGA * samtim + self.xref).sin_cos();

            x = ct * x1 - st * y1;
            y = st * x1 + ct * y1;

            let vcste1 = x - xsat;
            let vcste2 = y - ysat;
            let vcste3 = z - zsat;

            let vcses3 = self.b31 * vcste1 + self.b32 * vcste2 + self.b33 * vcste3;
            znorm = (vcste1.powi(2) + vcste2.powi(2) + vcste3.powi(2)).sqrt();

            let x3 = vcses3 / znorm;
            let umv = (x3.atan2((self.rfact - x3.powi(2)).sqrt())) - self.roasin;
            xlin = self.pic_center_line_num - umv / self.radlin;

            if i == 0 {
                samtim = self.time2;
                oldlin = xlin;
            }
        }

        let scnnum = ((oldlin + xlin) / 2.0 - 1.0) / self.numsen as f64;
        let scnfrc = (scnnum - self.scan1) / (self.scan2 - self.scan1);
        xlin = oldlin + scnfrc * (xlin - oldlin);
        samtim = self.time1 + self.tmpscl * (scnnum - self.scan1);
        let (new_xsat, new_ysat, new_zsat) = self.satvec(samtim);
        xsat = new_xsat;
        ysat = new_ysat;
        zsat = new_zsat;

        let cosa = x * xsat + y * ysat + z * zsat;
        let ctst = 0.0001 * RE * xht + RSQ;

        let mut xele = f64::NAN;

        if cosa >= ctst {
            let xsats1 = self.b11 * xsat + self.b12 * ysat + self.b13 * zsat;
            let ysats2 = self.b21 * xsat + self.b22 * ysat + self.b23 * zsat;
            let (st, ct) = (EMEGA * samtim + self.xref).sin_cos();
            x = ct * x1 - st * y1;
            y = st * x1 + ct * y1;
            let vcste1 = x - xsat;
            let vcste2 = y - ysat;
            let vcste3 = z - zsat;

            let vcses1 = self.b11 * vcste1 + self.b12 * vcste2 + self.b13 * vcste3;
            let vcses2 = self.b21 * vcste1 + self.b22 * vcste2 + self.b23 * vcste3;
            let vcses3 = self.b31 * vcste1 + self.b32 * vcste2 + self.b33 * vcste3;

            let xnorm = (znorm.powi(2) - vcses3.powi(2)).sqrt();
            let ynorm = (xsats1.powi(2) + ysats2.powi(2)).sqrt();
            znorm = (vcste1.powi(2) + vcste2.powi(2) + vcste3.powi(2)).sqrt();

            let x3 = vcses3 / znorm;
            let umv = x3.atan2((self.rfact - x3.powi(2)).sqrt()) - self.roasin;
            let (slin, clin) = umv.sin_cos();

            let u = self.rotm11 * clin + self.rotm13 * slin;
            let v = self.rotm21 * clin + self.rotm23 * slin;

            xele = self.picele
                + ((xsats1 * vcses2 - ysats2 * vcses1) / (xnorm * ynorm)).asin() / self.radele;
            xele += (v.atan2(u)) / self.radele;
            xele -= self.gamma + self.gamma_dot * samtim;
        }

        (xlin, xele)
    }
}
