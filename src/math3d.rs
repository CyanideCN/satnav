pub type Vec3 = [f64; 3];
pub type Mat3 = [[f64; 3]; 3];

pub fn dot(vector_a: &Vec3, vector_b: &Vec3) -> f64 {
    vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2]
}

pub fn magnitude(vector: &Vec3) -> f64 {
    dot(vector, vector).sqrt()
}

pub fn normalize(vector: &Vec3) -> Vec3 {
    let magnitude = magnitude(vector);
    if magnitude == 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [
            vector[0] / magnitude,
            vector[1] / magnitude,
            vector[2] / magnitude,
        ]
    }
}

pub fn cross(vector_a: &Vec3, vector_b: &Vec3) -> Vec3 {
    [
        vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1],
        vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2],
        vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0],
    ]
}

pub fn cross_normalized(vector_a: &Vec3, vector_b: &Vec3) -> Vec3 {
    normalize(&cross(vector_a, vector_b))
}

pub fn angle_between(vector_a: &Vec3, vector_b: &Vec3) -> f64 {
    let denominator = magnitude(vector_a) * magnitude(vector_b);
    if denominator == 0.0 {
        0.0
    } else {
        (dot(vector_a, vector_b) / denominator).acos()
    }
}

pub fn subtract(vector_a: &Vec3, vector_b: &Vec3) -> Vec3 {
    [
        vector_a[0] - vector_b[0],
        vector_a[1] - vector_b[1],
        vector_a[2] - vector_b[2],
    ]
}

pub fn add_scaled(vector_a: &Vec3, scale_a: f64, vector_b: &Vec3, scale_b: f64) -> Vec3 {
    [
        vector_a[0] * scale_a + vector_b[0] * scale_b,
        vector_a[1] * scale_a + vector_b[1] * scale_b,
        vector_a[2] * scale_a + vector_b[2] * scale_b,
    ]
}

pub fn mat3_mul_vec3(matrix: &Mat3, vector: &Vec3) -> Vec3 {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}

pub fn change_basis(basis_x: &Vec3, basis_y: &Vec3, basis_z: &Vec3, vector: &Vec3) -> Vec3 {
    [
        basis_x[0] * vector[0] + basis_y[0] * vector[1] + basis_z[0] * vector[2],
        basis_x[1] * vector[0] + basis_y[1] * vector[1] + basis_z[1] * vector[2],
        basis_x[2] * vector[0] + basis_y[2] * vector[1] + basis_z[2] * vector[2],
    ]
}
