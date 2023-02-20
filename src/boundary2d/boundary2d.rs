use anyhow::Result;
use std::collections::HashMap;
use nalgebra::base::*;
use imageproc::region_labelling::{connected_components, Connectivity};

pub struct Boundary2d {
    pub points: Vec<Vector2<f32>>,
    pub next_neighbor: Vec<usize>,
    pub prev_neighbor: Vec<usize>,
}

impl Boundary2d {
    pub fn from_naive_algorithm(mask: &image::GrayImage) -> Result<Boundary2d> {
        
        let mask_cc = connected_components(mask, Connectivity::Four, image::Luma([0u8]));

        let mut inds: HashMap<i32, usize> = HashMap::new();
        let mut edgs: Vec<[i32; 2]> = Vec::new();
        let mut nb_pts = 0;

        let mut points_wrk = Vec::new();
        let lin_len = (mask.width()+1) as i32;
        let cc_len = ((mask.width()+1) * (mask.height()+1)) as i32;

        for r_u32 in 0..(mask.height() + 1) {
            for c_u32 in 0..(mask.width() + 1) {
                // v00 | v10
                // v01 | v11

                let val00 = (if r_u32 != 0 && c_u32 != 0 {mask_cc.get_pixel(c_u32-1, r_u32-1).0[0]} else {0}) as i32;
                let val01 = (if r_u32 != 0 && c_u32 < mask_cc.width() {mask_cc.get_pixel(c_u32, r_u32-1).0[0]} else {0}) as i32;
                let val10 = (if r_u32 < mask_cc.height() && c_u32 != 0 {mask_cc.get_pixel(c_u32-1, r_u32).0[0]} else {0}) as i32;
                let val11 = (if r_u32 < mask_cc.height() && c_u32 < mask_cc.width() {mask_cc.get_pixel(c_u32, r_u32).0[0]} else {0}) as i32;
                
                let v00 = val00 != 0;
                let v01 = val01 != 0;
                let v10 = val10 != 0;
                let v11 = val11 != 0;

                if (val00 == val11 && !v01 && !v10) || (val10 == val01 && !v00 && !v11) {
                    return Err(anyhow::Error::msg("Adjacent corners not handled by naive boundary algorithm"));
                }

                if !((v00 && v01 && v10 && v11) || (!v00 && !v01 && !v10 && !v11)) {
                    let r_i32 = r_u32 as i32;
                    let c_i32 = c_u32 as i32;

                    let ind_cur = r_i32 * lin_len + c_i32;

                    points_wrk.push(Vector2::<f32>::new(c_i32 as f32, r_i32 as f32));
                    if val00 != 0 {
                        inds.insert(val00 * cc_len + ind_cur, nb_pts);
                        nb_pts = nb_pts+1;
                    }
                    if val01 != 0 && val01 != val00 {
                        inds.insert(val01 * cc_len + ind_cur, nb_pts);
                        nb_pts = nb_pts+1;
                    }
                    if val10 != 0 && val10 != val01 && val10 != val00 {
                        inds.insert(val10 * cc_len + ind_cur, nb_pts);
                        nb_pts = nb_pts+1;
                    }
                    if val11 != 0 && val11 != val10 && val11 != val01 && val11 != val00 {
                        inds.insert(val11 * cc_len + ind_cur, nb_pts);
                        nb_pts = nb_pts+1;
                    }

                    if !v00 && v01 {
                        let ind_up = (r_i32-1) * lin_len + c_i32;

                        let edg = [val01 * cc_len + ind_cur, val01 * cc_len + ind_up];

                        edgs.push(edg);
                    }

                    if !v01 && v11 {
                        let ind_right = r_i32 * lin_len + (c_i32+1);

                        let edg = [val11 * cc_len + ind_cur, val11 * cc_len + ind_right];

                        edgs.push(edg);
                    }

                    if !v11 && v10 {
                        let ind_down = (r_i32+1) * lin_len + c_i32;

                        let edg = [val10 * cc_len + ind_cur, val10 * cc_len + ind_down];

                        edgs.push(edg);
                    }

                    if !v10 && v00 {
                        let ind_left = r_i32 * lin_len + (c_i32-1);

                        let edg = [val00 * cc_len + ind_cur, val00 * cc_len + ind_left];

                        edgs.push(edg);
                    }
                }
            }
        }

        let mut next_neighbor_wrk = Vec::with_capacity(inds.len());
        let mut prev_neighbor_wrk = Vec::with_capacity(inds.len());

        next_neighbor_wrk.resize_with(inds.len(), Default::default);
        prev_neighbor_wrk.resize_with(inds.len(), Default::default);

        for [n1, n2] in edgs {
            let p1 = inds.get(&n1);
            let p2 = inds.get(&n2);

            match (p1, p2) {
                (Some(&i1), Some(&i2)) => {
                    next_neighbor_wrk[i1] = i2;
                    prev_neighbor_wrk[i2] = i1;
                },
                (None, _) => {return Err(anyhow::Error::msg("Index out of list"))},
                (_, None) => {return Err(anyhow::Error::msg("Index out of list"))},
            }
        }

        let bnd = Boundary2d{
            points: points_wrk,
            next_neighbor: next_neighbor_wrk,
            prev_neighbor: prev_neighbor_wrk
        };

        Ok(bnd)
    }
    
    pub fn from_marching_squares(mask: &image::GrayImage) -> Result<Vec<Boundary2d>> {
        
        let mask_cc = connected_components(mask, Connectivity::Eight, image::Luma([0u8]));
        
        struct MarchingData {
            points_march: HashMap<i32, Vector2<f32>>,
            edgs: Vec<[i32; 2]>,
        }
        
        let mut map_marching_data : HashMap<usize, MarchingData> = HashMap::new();

        let lin_len = 2 * (mask.width()+1) as i32;

        fn insert_and_connect(map_marching_data: &mut HashMap<usize, MarchingData>, 
                              lin_len: i32,
                              row: i32,
                              col: i32,
                              ind_cc: usize,
                              i0: usize,
                              i1: usize) -> () {
            let ind0 = (2 * row) * lin_len + (2 * col - 1);
            let ind1 = (2 * row - 1) * lin_len + (2 * col);
            let ind2 = (2 * row) * lin_len + (2 * col + 1);
            let ind3 = (2 * row + 1) * lin_len + (2 * col);
            let inds = [ind0, ind1, ind2, ind3];

            let pt0 = [col as f32 - 0.5, row as f32];
            let pt1 = [col as f32, row as f32 - 0.5];
            let pt2 = [col as f32 + 0.5, row as f32];
            let pt3 = [col as f32, row as f32 + 0.5];
            let pts = [pt0, pt1, pt2, pt3];
            
            match map_marching_data.get_mut(&ind_cc) {
                Some(map_cc) => {
                    map_cc.points_march.insert(inds[i0], Vector2::<f32>::new(pts[i0][0], pts[i0][1]));
                    map_cc.points_march.insert(inds[i1], Vector2::<f32>::new(pts[i1][0], pts[i1][1]));
                    map_cc.edgs.push([inds[i0], inds[i1]]);
                }
                None => {
                    let mut map_cc = MarchingData {
                        points_march : HashMap::new(),
                        edgs : Vec::new(),
                    };
                    map_cc.points_march.insert(inds[i0], Vector2::<f32>::new(pts[i0][0], pts[i0][1]));
                    map_cc.points_march.insert(inds[i1], Vector2::<f32>::new(pts[i1][0], pts[i1][1]));
                    map_cc.edgs.push([inds[i0], inds[i1]]);
                    map_marching_data.insert(ind_cc, map_cc);
                }
            }
        }

        for r_u32 in 0..(mask.height() + 1) {
            for c_u32 in 0..(mask.width() + 1) {
                // v00 1 v10
                //  0  +  2
                // v01 3 v11

                let val00 = (if r_u32 != 0 && c_u32 != 0 {mask_cc.get_pixel(c_u32-1, r_u32-1).0[0]} else {0}) as usize;
                let val01 = (if r_u32 < mask_cc.height() && c_u32 != 0 {mask_cc.get_pixel(c_u32-1, r_u32).0[0]} else {0}) as usize;
                let val10 = (if r_u32 != 0 && c_u32 < mask_cc.width() {mask_cc.get_pixel(c_u32, r_u32-1).0[0]} else {0}) as usize;
                let val11 = (if r_u32 < mask_cc.height() && c_u32 < mask_cc.width() {mask_cc.get_pixel(c_u32, r_u32).0[0]} else {0}) as usize;
                
                let v00 = val00 != 0;
                let v01 = val01 != 0;
                let v10 = val10 != 0;
                let v11 = val11 != 0;

                let r_i32 = r_u32 as i32;
                let c_i32 = c_u32 as i32;

                // 1 corner: 4 cases
                if v00 && !v10 && !v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val00, 1, 0);
                }
                if !v00 && v10 && !v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val10, 2, 1);
                }
                if !v00 && !v10 && v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val11, 3, 2);
                }
                if !v00 && !v10 && !v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val01, 0, 3);
                }
                                                                                  
                // 2 consecutives : 4 cases                                          
                if v00 && v10 && !v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val00, 2, 0);
                }
                if !v00 && v10 && v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val10, 3, 1);
                }
                if !v00 && !v10 && v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val11, 0, 2);
                }
                if v00 && !v10 && !v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val01, 1, 3);
                }
                
                // 3 consecutives : 4 cases                                          
                if v00 && v10 && v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val00, 3, 0);
                }
                if !v00 && v10 && v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val10, 0, 1);
                }
                if v00 && !v10 && v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val11, 1, 2);
                }
                if v00 && v10 && !v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val01, 2, 3);
                }
                
                // 2 opposite (ambiguous case) : 2 cases                                        
                if v00 && !v10 && v11 && !v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val00, 3, 0);
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val00, 1, 2);
                }                                                                    
                if !v00 && v10 && !v11 && v01 {
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val10, 0, 1);
                    insert_and_connect(&mut map_marching_data, lin_len, r_i32, c_i32, val10, 2, 3);
                }
            }
        }
        
        let mut vec_bnd = Vec::new();
        for (_, marching_data) in map_marching_data {
            let mut inds = HashMap::new();
            let mut points_wrk = Vec::new();
            marching_data
                .points_march
                .iter()
                .enumerate()
                .fold((&mut inds, &mut points_wrk), 
                      |(inds, points_wrk), (cpt, (num, &pt))| {
                          inds.insert(num, cpt); 
                          points_wrk.push(pt); 
                          (inds, points_wrk)
                      });

            let mut next_neighbor_wrk = Vec::with_capacity(inds.len());
            let mut prev_neighbor_wrk = Vec::with_capacity(inds.len());

            next_neighbor_wrk.resize_with(inds.len(), Default::default);
            prev_neighbor_wrk.resize_with(inds.len(), Default::default);

            for [n1, n2] in marching_data.edgs {
                let p1 = inds.get(&n1);
                let p2 = inds.get(&n2);

                match (p1, p2) {
                    (Some(&i1), Some(&i2)) => {
                        next_neighbor_wrk[i1] = i2;
                        prev_neighbor_wrk[i2] = i1;
                    },
                    (None, _) => {return Err(anyhow::Error::msg("Index out of list"))},
                    (_, None) => {return Err(anyhow::Error::msg("Index out of list"))},
                }
            }
            
            let bnd = Boundary2d{
                points: points_wrk,
                next_neighbor: next_neighbor_wrk,
                prev_neighbor: prev_neighbor_wrk
            };
            vec_bnd.push(bnd);
        }
        
        Ok(vec_bnd)
    }

    fn get_loops_ind(&self) -> Vec<Vec<usize>> {
        let mut used = vec![false; self.points.len()];
        
        let mut loops_ind = Vec::<Vec<usize>>::new();

        for i in 0..used.len() {
            if !used[i] {
                let mut loop_ind = Vec::new();

                let mut ind_cur = i as usize;

                while !used[ind_cur] {
                    loop_ind.push(ind_cur);
                    used[ind_cur] = true;
                    ind_cur = self.next_neighbor[ind_cur];
                }

                loops_ind.push(loop_ind);
            }
        }

        loops_ind
    }

    pub fn get_loops_pts(&self) -> Vec<Vec<Vector2<f32>>> {
        self.get_loops_ind()
            .into_iter()
            .map(
                |v| 
                v.into_iter()
                .map(|ind| self.points[ind]).collect()).collect()
    }
}


