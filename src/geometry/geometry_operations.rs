use nalgebra::base::*;

pub fn circle_center(
    pt1: &Vector2<f32>,
    pt2: &Vector2<f32>,
    pt3: &Vector2<f32>,
) -> Option<Vector2<f32>> {
    // I'm sure there is a justified equation behind those lines, but I don't remember it...

    let mat = Matrix3x2::new(
        pt2[0] - pt1[0],
        pt2[1] - pt1[1],
        pt3[0] - pt2[0],
        pt3[1] - pt2[1],
        pt1[0] - pt3[0],
        pt1[1] - pt3[1],
    );

    let b = Vector3::new(
        0.5 * (pt2 - pt1).norm_squared() + (pt2 - pt1).dot(pt1),
        0.5 * (pt3 - pt2).norm_squared() + (pt3 - pt2).dot(pt2),
        0.5 * (pt1 - pt3).norm_squared() + (pt1 - pt3).dot(pt3),
    );

    let mat_mod = mat.transpose() * mat;
    let b_mod = mat.transpose() * b;

    mat_mod.lu().solve(&b_mod)
}

/// Computes sphere center associated to the 4 given points
pub fn sphere_center(
    pt1: &Vector3<f32>,
    pt2: &Vector3<f32>,
    pt3: &Vector3<f32>,
    pt4: &Vector3<f32>,
) -> Option<Vector3<f32>> {
    let vec_0_1 = pt2 - pt1;
    let vec_1_2 = pt3 - pt2;
    let vec_2_0 = pt1 - pt3;
    let vec_3_0 = pt1 - pt4;
    let vec_3_1 = pt2 - pt4;
    let vec_3_2 = pt3 - pt4;

    #[rustfmt::skip]
    let mat_slv = Matrix6x3::new(
        vec_0_1[0], vec_0_1[1], vec_0_1[2], 
        vec_1_2[0], vec_1_2[1], vec_1_2[2], 
        vec_2_0[0], vec_2_0[1], vec_2_0[2], 
        vec_3_0[0], vec_3_0[1], vec_3_0[2], 
        vec_3_1[0], vec_3_1[1], vec_3_1[2], 
        vec_3_2[0], vec_3_2[1], vec_3_2[2], 
    );

    let sqn0 = pt1.norm_squared();
    let sqn1 = pt2.norm_squared();
    let sqn2 = pt3.norm_squared();
    let sqn3 = pt4.norm_squared();

    let vec_slv = Matrix6x1::new(
        0.5 * (sqn1 - sqn0),
        0.5 * (sqn2 - sqn1),
        0.5 * (sqn0 - sqn2),
        0.5 * (sqn0 - sqn3),
        0.5 * (sqn1 - sqn3),
        0.5 * (sqn2 - sqn3),
    );
    let mat_slv_mod = mat_slv.transpose() * mat_slv;
    let vec_slv_mod = mat_slv.transpose() * vec_slv;

    mat_slv_mod.lu().solve(&vec_slv_mod)
}
