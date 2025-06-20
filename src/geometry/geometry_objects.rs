use crate::geometry::geometry_operations;
use anyhow::Result;
use nalgebra::base::*;

#[derive(Copy, Clone)]
pub struct Circle {
    pub center: Vector2<f32>,
    pub radius: f32,
}

impl Circle {
    pub fn from(pt1: &Vector2<f32>, pt2: &Vector2<f32>, pt3: &Vector2<f32>) -> Result<Circle> {
        let ctr = geometry_operations::circle_center(pt1, pt2, pt3)
            .ok_or(anyhow::Error::msg("Could not compute circle center"))?;
        let rad = (ctr - pt1).norm();

        Ok(Circle {
            center: ctr,
            radius: rad,
        })
    }
}

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Vector3<f32>,
    pub radius: f32,
}

impl Sphere {
    pub fn from(
        pt1: &Vector3<f32>,
        pt2: &Vector3<f32>,
        pt3: &Vector3<f32>,
        pt4: &Vector3<f32>,
    ) -> Result<Sphere> {
        let ctr = geometry_operations::sphere_center(pt1, pt2, pt3, pt4)
            .ok_or(anyhow::Error::msg("Could not compute sphere center"))?;
        let rad = (ctr - pt1).norm();

        Ok(Sphere {
            center: ctr,
            radius: rad,
        })
    }
}
