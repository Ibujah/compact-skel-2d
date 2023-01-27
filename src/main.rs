use clap::Parser;
use anyhow::Result;
use std::time::Instant;

use image::GenericImageView;
use imageproc;

use compact_skel_2d::skeleton2d;
use compact_skel_2d::boundary2d;

#[derive(Parser)]
struct Cli {
    #[arg(default_value="./ressources/rat.png", long = "imgfile")]
    image_path: std::path::PathBuf,

    #[arg(default_value="1.0", long = "epsilon")]
    epsilon: f32,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    
    let image_path_str = args.image_path.to_str().unwrap_or("");
    let epsilon = args.epsilon;
    
    let mask = image::open(image_path_str)?;
    let mask_luma8 = mask.as_luma8().ok_or(anyhow::Error::msg("Could not convert to luma8"))?;
    let mask_rgb = 
        image::RgbImage::from_vec(
            mask.width(), 
            mask.height(), 
            mask.as_bytes()
            .to_vec()
            .into_iter()
            .fold(
                Vec::<u8>::new(),
                |mut vec, x| 
                {
                    let v = if x >= 125 {255} else {0};
                    vec.push(v);vec.push(v);vec.push(v);
                    vec
                }
                )
            ).ok_or(anyhow::Error::msg("Could not convert gray to rgb"))?;
    let mut mask_rgb = image::DynamicImage::ImageRgb8(mask_rgb);

    
    let now = Instant::now();
    println!("Boundary computation");
    let vec_bnd = boundary2d::Boundary2d::from_marching_squares(mask_luma8)?;
    for bnd in vec_bnd {
        println!("Skeleton initialization");
        let mut skel = skeleton2d::Skeleton::init(&bnd.points, &bnd.next_neighbor, &bnd.prev_neighbor);
        println!("Delaunay triangles computation");
        skel.compute_delaunay_triangles()?;
        println!("Finding first triangle");
        let ind_tri = skeleton2d::find_first_in(&mut skel)?;
        println!("Finding nearest maximum");
        let ind_max = skeleton2d::find_nearest_maximum(&mut skel, ind_tri)?;
        println!("Skeleton propagation");
        skeleton2d::propagate_from(&mut skel, ind_max, epsilon)?;
        let duration = now.elapsed();

        println!("Skeleton computed in {}ms", duration.as_millis());

        let circle = skel.get_circle(ind_max)?;
        imageproc::drawing::draw_hollow_circle_mut(&mut mask_rgb, 
                                                   ((circle.center[0]+0.5) as i32, (circle.center[1]+0.5) as i32), 
                                                   (circle.radius+0.5) as i32, 
                                                   image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]));

        let edges_ind = skel.get_edges_ind();

        for [i1, i2] in edges_ind {
            let cir1 = skel.get_circle(i1)?;
            let cir2 = skel.get_circle(i2)?;
            let ctr1 = cir1.center;
            let ctr2 = cir2.center;
            imageproc::drawing::draw_line_segment_mut(&mut mask_rgb, 
                                                      (ctr1[0], ctr1[1]), 
                                                      (ctr2[0], ctr2[1]), 
                                                      image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]));
        }
    }

    image::save_buffer("exp.png", mask_rgb.as_bytes(), mask_rgb.width(), mask_rgb.height(), mask_rgb.color())?;

    Ok(())
}
