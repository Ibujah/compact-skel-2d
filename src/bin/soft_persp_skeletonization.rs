use anyhow::Result;
use clap::Parser;
use env_logger;
use log;
use std::time::Instant;

use image::GenericImageView;
use imageproc;
use nalgebra::base::*;

use compact_skel_2d::boundary2d::boundary2d;
use compact_skel_2d::skeleton2d::skeleton_persp::{self, SkeletonPersp};

#[derive(Parser)]
struct Cli {
    #[arg(default_value = "./ressources/rat.png", long = "imgin")]
    image_in_path: std::path::PathBuf,

    #[arg(default_value = "./ressources/rat_skel.png", long = "imgout")]
    image_out_path: std::path::PathBuf,

    #[arg(default_value = "800", long = "fu")]
    fu: f32,

    #[arg(default_value = "800", long = "fv")]
    fv: f32,

    #[arg(default_value = "400", long = "cu")]
    cu: f32,

    #[arg(default_value = "300", long = "cv")]
    cv: f32,

    #[arg(default_value = "1.0", long = "epsilon")]
    epsilon: f32,

    #[arg(default_value = "./ressources/nod.txt", long = "nodout")]
    nod_out_path: std::path::PathBuf,

    #[arg(default_value = "./ressources/edg.txt", long = "edgout")]
    edg_out_path: std::path::PathBuf,

    #[arg(default_value = "./ressources/del.txt", long = "delout")]
    del_out_path: std::path::PathBuf,

    #[arg(default_value = "./ressources/bnd.txt", long = "bndout")]
    bnd_out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let image_in_path_str = args.image_in_path.to_str().unwrap_or("");
    let image_out_path_str = args.image_out_path.to_str().unwrap_or("");
    let nod_out_path_str = args.nod_out_path.to_str().unwrap_or("");
    let edg_out_path_str = args.edg_out_path.to_str().unwrap_or("");
    let del_out_path_str = args.del_out_path.to_str().unwrap_or("");
    let bnd_out_path_str = args.bnd_out_path.to_str().unwrap_or("");
    let epsilon = args.epsilon;
    let fu = args.fu;
    let fv = args.fv;
    let cu = args.cu;
    let cv = args.cv;

    let mask = image::open(image_in_path_str)?;
    let mask_luma8 = mask
        .as_luma8()
        .ok_or(anyhow::Error::msg("Could not convert to luma8"))?;
    let mask_rgb = image::RgbImage::from_vec(
        mask.width(),
        mask.height(),
        mask.as_bytes()
            .to_vec()
            .into_iter()
            .fold(Vec::<u8>::new(), |mut vec, x| {
                let v = if x >= 125 { 255 } else { 0 };
                vec.push(v);
                vec.push(v);
                vec.push(v);
                vec
            }),
    )
    .ok_or(anyhow::Error::msg("Could not convert gray to rgb"))?;
    let mut mask_rgb = image::DynamicImage::ImageRgb8(mask_rgb);

    let now = Instant::now();
    log::info!("Boundary computation");
    let vec_bnd = boundary2d::Boundary2d::from_marching_squares(mask_luma8)?;
    let mut skel = SkeletonPersp::new(fu, fv, cu, cv);
    for bnd in vec_bnd {
        log::info!("Add boundary to skeleton and find first triangle");
        let ind_tri = skeleton_persp::append_and_find_first(
            &mut skel,
            &bnd.points,
            &bnd.next_neighbor,
            &bnd.prev_neighbor,
        )?;
        log::info!("Finding nearest maximum");
        let ind_max = skeleton_persp::find_nearest_maximum(&mut skel, ind_tri)?;
        log::info!("Skeleton propagation");
        skeleton_persp::propagate_from(&mut skel, ind_max, epsilon)?;
    }
    let duration = now.elapsed();
    log::info!("Skeleton computed in {}ms", duration.as_millis());

    let edges_ind = skel.get_edges_ind();

    for [i1, i2] in edges_ind {
        let sph1 = skel.get_sphere(i1)?;
        let sph2 = skel.get_sphere(i2)?;
        let ctr1 = sph1.center;
        let ctr2 = sph2.center;
        let proj_ctr1 = Vector2::new(fu * ctr1[0] / ctr1[2] + cu, fv * ctr1[1] / ctr1[2] + cv);
        let proj_ctr2 = Vector2::new(fu * ctr2[0] / ctr2[2] + cu, fv * ctr2[1] / ctr2[2] + cv);

        imageproc::drawing::draw_line_segment_mut(
            &mut mask_rgb,
            (proj_ctr1[0], proj_ctr1[1]),
            (proj_ctr2[0], proj_ctr2[1]),
            image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
        );
    }

    image::save_buffer(
        image_out_path_str,
        mask_rgb.as_bytes(),
        mask_rgb.width(),
        mask_rgb.height(),
        mask_rgb.color(),
    )?;

    skeleton_persp::export_skel_txt(
        &skel,
        nod_out_path_str,
        edg_out_path_str,
        del_out_path_str,
        bnd_out_path_str,
    )?;

    Ok(())
}
