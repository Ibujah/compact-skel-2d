use nalgebra::base::*;

pub fn circle_center(pt1: &Vector2<f32>, pt2: &Vector2<f32>, pt3: &Vector2<f32>) -> Option<Vector2<f32>> {
    // I'm sure there is a justified equation behind those lines, but I don't remember it...

    let mat = Matrix3x2::new(
        pt2[0] - pt1[0], pt2[1] - pt1[1], 
        pt3[0] - pt2[0], pt3[1] - pt2[1], 
        pt1[0] - pt3[0], pt1[1] - pt3[1], 
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
