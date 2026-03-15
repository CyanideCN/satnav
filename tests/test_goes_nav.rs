#[cfg(test)]
mod goes_tests {
    use satnav::goes_nav::GOESNavigation;

    const NAV_BLK: [i32; 128] = [
        1196377427, 3292276, 213100, 1, 921002, 0, 4216558, 340, 141, 194546, 324369, 99787,
        894843, 104652, 7285, 599906, 200000, 801821, 182230, 15288, 2537, -2900, 0, 0, 0, 0, 0,
        453359, 0, 0, 154, 213104, 60, 3199410, 1751, 214702, 66, 3129676, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    #[test]
    fn test_calc_subpoint() {
        let nav = GOESNavigation::new(&NAV_BLK);
        let (sub_lat, sub_lon) = nav.get_subpoint();
        assert!((sub_lat - 0.11960715927181872).abs() < 1e-7);
        assert!((sub_lon - -112.24144895816976).abs() < 1e-7);
    }

    #[test]
    fn test_lineele_lonlat_conv() {
        let nav = GOESNavigation::new(&NAV_BLK);
        let line = 5000.;
        let ele = 5000.;
        let (lat, lon) = nav.image_to_geographic(line, ele);
        let (nlin, nele) = nav.geographic_to_image(lat, lon);
        assert!((nlin - line).abs() < 0.5);
        assert!((nele - ele).abs() < 0.5);
    }
}
