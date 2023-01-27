use anyhow::Result;
use nalgebra::base::*;
use crate::geometry::geometry_operations;

#[derive(Copy, Clone)]
pub struct Circle{
    pub center: Vector2<f32>,
    pub radius: f32,
}

impl Circle{
    pub fn from(pt1: &Vector2<f32>, pt2: &Vector2<f32>, pt3: &Vector2<f32>) -> Result<Circle> {
        let ctr = geometry_operations::circle_center(pt1, pt2, pt3).ok_or(anyhow::Error::msg("Could not compute circle center"))?;
        let rad = (ctr - pt1).norm();

        Ok(Circle{
            center: ctr,
            radius: rad,
        })
    }
}
