use anyhow::Result;
use std::collections::HashMap;
use delaunator::{Point, triangulate};
use nalgebra::base::*;

use crate::geometry::geometry_objects::Circle;

pub struct Skeleton{
    pub boundary_points: Vec<Vector2<f32>>,
    pub next_neighbor: Vec<usize>,
    pub prev_neighbor: Vec<usize>,

    pub delaunay_triangles: Vec<[usize; 3]>,
    pub delaunay_neighbors: HashMap<usize, [Option<usize>; 3]>,
    pub skel_circles: HashMap<usize, Circle>,

    pub final_skeleton: Vec<bool>,
}

impl Skeleton {
    pub fn init(boundary_points: &Vec<Vector2<f32>>, next_neighbor: &Vec<usize>, prev_neighbor: &Vec<usize>) -> Skeleton {
        Skeleton{
            boundary_points: boundary_points.clone(), 
            next_neighbor: next_neighbor.clone(), 
            prev_neighbor: prev_neighbor.clone(), 
            delaunay_triangles: Vec::new(),
            delaunay_neighbors : HashMap::new(),
            skel_circles: HashMap::new(),
            final_skeleton: Vec::new(),
        }
    }

    pub fn compute_delaunay_triangles(&mut self) -> Result<()> {

        let points: Vec<Point> = 
            self
            .boundary_points
            .iter()
            .map(|v| 
                 Point{x: v[0] as f64, y: v[1] as f64}).collect();

        let result = triangulate(&points);

        let mut iter = result.triangles.iter();

        loop {
            let v1 = iter.next();
            let v2 = iter.next();
            let v3 = iter.next();

            match (v1, v2, v3) {
                (Some(&i1), Some(&i2), Some(&i3)) => {
                    self.delaunay_triangles.push([i1, i2, i3]);
                },
                (None, None, None) => 
                    break,
                (_, _, _) => 
                    return Err(anyhow::Error::msg("Error collecting delaunay triangles")),
            }
        }

        self.final_skeleton = vec![false; self.delaunay_triangles.len()];
    
        Ok(())
    }
    
    fn find_neighbors(&self, ind_triangle: usize) -> [Option<usize>; 3] {
        let [i1, i2, i3] = self.delaunay_triangles[ind_triangle];
        
        let nei = 
            self.delaunay_triangles
            .iter()
            .enumerate()
            .fold([None, None, None],
            |[n1, n2, n3], (ind, v)|
            {
                let i1_in = v[0] == i1 || v[1] == i1 || v[2] == i1;
                let i2_in = v[0] == i2 || v[1] == i2 || v[2] == i2;
                let i3_in = v[0] == i3 || v[1] == i3 || v[2] == i3;

                let new_nei = 
                    if i1_in && i2_in && !i3_in {
                        [n1, n2, Some(ind)]
                    }
                    else if i1_in && !i2_in && i3_in {
                        [n1, Some(ind), n3]
                    }
                    else if !i1_in && i2_in && i3_in {
                        [Some(ind), n2, n3]
                    }
                    else {
                        [n1, n2, n3]
                    };

                new_nei
            }
        );

        nei
    }

    pub fn get_neighbors(&mut self, ind_triangle: usize) -> [Option<usize>; 3] {
        match self.delaunay_neighbors.get(&ind_triangle) {
            Some(val) => *val,
            None => {
                let nei = self.find_neighbors(ind_triangle);
                self.delaunay_neighbors.insert(ind_triangle, nei);
                nei
            }
        }
    }

    pub fn get_non_crossing_neighbors(&mut self, ind_triangle: usize) -> [Option<usize>; 3] {
        let all_nei = self.get_neighbors(ind_triangle);
        let ind_tri = self.delaunay_triangles[ind_triangle];
        let mut kept_nei = [None, None, None];
        for i in 0..3 {
            kept_nei[i] = 
                match all_nei[i] {
                    Some(val) => {
                        let i1 = ind_tri[(i+1)%3];
                        let i2 = ind_tri[(i+2)%3];

                        if self.next_neighbor[i2] == i1 {
                            None
                        }
                        else {
                            Some(val)
                        }
                    },
                    None => {
                        None
                    }
                }
        }

        kept_nei
    }

    pub fn get_non_covered_neighbors(&mut self, ind_triangle: usize, epsilon: f32) -> Result<[Option<usize>; 3]> {
        let all_nei = self.get_neighbors(ind_triangle);
        let ind_tri = self.delaunay_triangles[ind_triangle];
        let cir = self.get_circle(ind_triangle)?;
        let mut kept_nei = [None, None, None];
        for i in 0..3 {
            kept_nei[i] = 
                match all_nei[i] {
                    Some(val) => {
                        let i1 = ind_tri[(i+1)%3];
                        let i2 = ind_tri[(i+2)%3];
                        let mut icur = i2;
                        let res;
                        loop {
                            if self.next_neighbor[icur] == i1 {
                                res = None;
                                break
                            }
                            else {
                                let pt = self.boundary_points[icur];
                                let dist = (pt - cir.center).norm();
                                if dist < cir.radius + epsilon {
                                    icur = self.next_neighbor[icur];
                                }
                                else {
                                    res = Some(val);
                                    break;
                                }
                                if icur == i2 {
                                    res = Some(val);
                                    break;
                                }
                            } 
                        }
                        res
                    },
                    None => {
                        None
                    }
                }
        }

        Ok(kept_nei)
    }

    pub fn get_circle(&mut self, ind_triangle: usize) -> Result<Circle> {
        match self.skel_circles.get(&ind_triangle) {
            Some(val) => Ok(*val),
            None => {
                let [ind1, ind2, ind3] = self.delaunay_triangles[ind_triangle];

                let pt1 = self.boundary_points[ind1];
                let pt2 = self.boundary_points[ind2];
                let pt3 = self.boundary_points[ind3];

                let res = Circle::from(&pt1, &pt2, &pt3)?;

                self.skel_circles.insert(ind_triangle, res);
                Ok(res)
            }
        }
    }

    pub fn get_edges_ind(&mut self) -> Vec<[usize; 2]> {
        let mut edges_ind = Vec::new();

        for i in 0..self.final_skeleton.len() {
            if self.final_skeleton[i] {
                let nei = self.get_non_crossing_neighbors(i);
                for j in 0..3 {
                    match nei[j] {
                        Some(v) => {
                            if v > i && self.final_skeleton[v] {
                                edges_ind.push([i, v]);
                            }
                        }
                        None => ()
                    }
                }
            }
        }
        edges_ind
    }
}

pub fn find_first_in(skeleton: &mut Skeleton) -> Result<usize> {
    let corner_box = 
        skeleton.boundary_points
        .iter()
        .fold(None, 
              |corner: Option<Vector2<f32>>, v|
              {
                  match corner {
                      Some(cor) =>
                      {
                          let min_x = if cor[0] < v[0] {cor[0]} else {v[0]};
                          let min_y = if cor[1] < v[1] {cor[1]} else {v[1]};
                          Some(Vector2::<f32>::new(min_x, min_y))
                      },
                      None => {Some(*v)}
                  }
              }
             );
    
    let corner = corner_box.ok_or(anyhow::Error::msg("Empty boundary"))?;

    let (ind_clo, _) = 
        skeleton.boundary_points
        .iter()
        .enumerate()
        .fold((0, 0.0),
        |(ind_min, dist_min), (ind_cur, v)|
        {
            let dist_cur = (v - corner).norm();
            if ind_cur == 0 || dist_cur < dist_min {
                (ind_cur, dist_cur)
            }
            else {
                (ind_min, dist_min)
            }
        }
        );

    let ind_next = skeleton.next_neighbor[ind_clo];
    
    let ind_tri = 
        skeleton.delaunay_triangles
        .iter()
        .position(
            |&[n1, n2, n3]|
            {
                let config1 = n2 == ind_clo && n1 == ind_next;
                let config2 = n3 == ind_clo && n2 == ind_next;
                let config3 = n1 == ind_clo && n3 == ind_next;

                config1 || config2 || config3
            }
            ).ok_or(anyhow::Error::msg("Could not find first triangle"))?;

    Ok(ind_tri)
}

pub fn find_nearest_maximum(skeleton: &mut Skeleton, ind_triangle: usize) -> Result<usize> {
    let mut ind = ind_triangle;
    let mut cir = skeleton.get_circle(ind_triangle)?;

    loop {
        let nei = skeleton.get_non_crossing_neighbors(ind);
        let upd: bool;
        
        (upd, ind, cir) = 
            nei.iter()
            .fold(Ok((false, ind, cir)),
            |cur:Result<(bool, usize, Circle)>, n| {
                let (upd, cur_ind, cur_cir) = cur?;
                match n {
                    Some(i) => {
                        let cir = skeleton.get_circle(*i)?;
                        if cir.radius > cur_cir.radius {
                            Ok((true, *i, cir))
                        }
                        else {
                            Ok((upd, cur_ind, cur_cir))
                        }
                    },
                    None => {
                        Ok((upd, cur_ind, cur_cir))
                    }
                }
            })?;

        if !upd {
            break
        }
    }
    Ok(ind)
}

pub fn propagate_from(skeleton: &mut Skeleton, ind_triangle: usize, epsilon: f32) -> Result<()> {
    let mut to_propagate = vec![ind_triangle];

    loop {
        match to_propagate.pop() {
            Some(ind) => {
                if !skeleton.final_skeleton[ind] {
                    let nei = skeleton.get_non_covered_neighbors(ind, epsilon)?;
                    for n in nei {
                        match n {
                            Some(v) => to_propagate.push(v),
                            None => ()
                        }
                    }
                    skeleton.final_skeleton[ind] = true;
                }
            }
            None => {
                break
            }
        }
    }

    Ok(())
}





