use clap::Parser;
use image::{Rgb, RgbImage};
use num_complex::Complex64;
use std::thread;

pub mod polynomial;
use crate::polynomial::Polynome;

const WIDTH: u32 = 1000;
const HEIGHT: u32 = 1000;
const MAX_ITER: u32 = 60;
const ROOT_MERGE_DIST: f64 = 1e-4;

const X_MIN: f64 = -2.0;
const X_MAX: f64 = 2.0;
const Y_MIN: f64 = -2.0;
const Y_MAX: f64 = 2.0;

#[derive(Parser)]
struct Args {
    /// Coefficients du polynôme, du degré le plus bas au plus haut
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true)]
    coeffs: Option<Vec<f64>>,

    /// Racines du polynôme sous forme de paires re,im : re1,im1,re2,im2,...
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true)]
    roots: Option<Vec<f64>>,

    #[arg(long, short, default_value = "wada.png")]
    path: String,
}

struct ThreadResult {
    row_start: u32,
    iterations: Vec<u32>,
    local_root_idx: Vec<Option<usize>>, // None = non convergé
    local_roots: Vec<Complex64>,
}

enum RootResolution {
    Known,
    Discovered { mapping: Vec<Vec<usize>> },
}

impl RootResolution {
    fn global_index(&self, thread_idx: usize, local_idx: usize) -> usize {
        match self {
            RootResolution::Known => local_idx,
            RootResolution::Discovered { mapping } => mapping[thread_idx][local_idx],
        }
    }
}
///Regroupe les racines trouvées par chaque thread, fusionne les racines trop proches, seuil de
///tolérance : merge_dist
pub fn cluster_roots(
    per_thread_roots: &[Vec<Complex64>],
    merge_dist: f64,
) -> (Vec<Complex64>, Vec<Vec<usize>>) {
    let mut global_roots: Vec<Complex64> = Vec::new();
    let mut mapping: Vec<Vec<usize>> = Vec::with_capacity(per_thread_roots.len());

    for thread_roots in per_thread_roots {
        let mut thread_mapping = Vec::with_capacity(thread_roots.len());
        for &r in thread_roots {
            let existing = global_roots
                .iter()
                .position(|g| (r - g).norm() < merge_dist);
            let global_idx = match existing {
                Some(idx) => idx,
                None => {
                    global_roots.push(r);
                    global_roots.len() - 1
                }
            };
            thread_mapping.push(global_idx);
        }
        mapping.push(thread_mapping);
    }

    (global_roots, mapping)
}

///Renvoie une valeur rgb en fonction pour une racine d'un polynome de degré deg
fn color_for(root_idx: usize, iterations: u32, deg: usize) -> Rgb<u8> {
    let palette = [
        (255u8, 80u8, 80u8),   // rouge
        (255u8, 160u8, 80u8),  // orange
        (255u8, 220u8, 80u8),  // jaune
        (180u8, 255u8, 80u8),  // vert-jaune
        (80u8, 255u8, 120u8),  // vert
        (80u8, 255u8, 200u8),  // vert-cyan
        (80u8, 230u8, 230u8),  // cyan
        (80u8, 180u8, 255u8),  // bleu clair
        (100u8, 140u8, 255u8), // bleu
        (140u8, 100u8, 255u8), // indigo
        (200u8, 100u8, 255u8), // violet
        (255u8, 100u8, 220u8), // magenta
        (255u8, 100u8, 160u8), // rose
        (200u8, 200u8, 200u8), // gris clair
        (255u8, 255u8, 255u8), // blanc
        (150u8, 150u8, 80u8),  // olive
    ];
    let (r, g, b) = palette[(root_idx % deg) % palette.len()];
    let shade = 1.0 - (iterations as f32 / MAX_ITER as f32) * 0.7;
    Rgb([
        (r as f32 * shade) as u8,
        (g as f32 * shade) as u8,
        (b as f32 * shade) as u8,
    ])
}

fn parse_complex_pairs(flat: &[f64]) -> Vec<Complex64> {
    flat.chunks(2)
        .map(|pair| Complex64::new(pair[0], pair.get(1).copied().unwrap_or(0.0)))
        .collect()
}

/// Génère une image PNG de la fractale de Newton associée à un polynôme.
///
/// Le polynôme peut être fourni de deux façons (mutuellement exclusives) :
///
/// - `--coeffs re0,im0,re1,im1,...` : coefficients du degré le plus bas au plus haut
/// - `--roots re0,im0,re1,im1,...` : racines connues du polynôme
///
/// Le chemin de sortie peut être précisé avec `--path` (par défaut `wada.png`).
///
/// # Examples
///
/// ```bash
/// # z^3 - 1 via ses coefficients
/// cargo run --release -- --coeffs=-1,0,0,0,0,0,1,0
/// wada --coeffs -1,0,0,0,0,0,1,0
///
/// # z^3 - 1 via ses racines directement
/// cargo run --release -- --roots 1,0,-0.5,0.87,-0.5,-0.87 --path racines.png
/// wada -p racines.png --roots 1,0,-0.5,0.87,-0.5,-0.87
/// ```
fn main() {
    let args = Args::parse();
    let path = args.path;

    let known_roots: Option<Vec<Complex64>> = args.roots.as_ref().map(|r| parse_complex_pairs(r));

    let poly = if let Some(roots) = &known_roots {
        if roots.is_empty() {
            eprintln!("Il faut au moins une racine.");
            std::process::exit(1);
        }
        Polynome::from_roots(roots)
    } else {
        let coeffs_flat = args
            .coeffs
            .unwrap_or_else(|| vec![-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);

        if coeffs_flat.len() % 2 != 0 {
            eprintln!("--coeffs doit contenir un nombre pair de valeurs (paires re,im).");
            std::process::exit(1);
        }

        let coeffs = parse_complex_pairs(&coeffs_flat);

        if coeffs.len() < 2 {
            eprintln!("Il faut au moins un polynôme de degré 1 (2 coefficients).");
            std::process::exit(1);
        }
        Polynome::new(coeffs)
    };

    let dpoly = poly.derivative();

    let n_threads = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let rows_per_thread = HEIGHT.div_ceil(n_threads as u32);

    let mut raw_results: Vec<ThreadResult> = Vec::new();

    thread::scope(|s| {
        let mut handles = Vec::new();

        for t in 0..n_threads {
            let row_start = t as u32 * rows_per_thread;
            let row_end = (row_start + rows_per_thread).min(HEIGHT);
            if row_start >= row_end {
                continue;
            }

            let poly_ref = &poly;
            let dpoly_ref = &dpoly;
            let known_roots_ref = known_roots.as_deref();

            let handle = s.spawn(move || {
                let pixel_count = ((row_end - row_start) * WIDTH) as usize;
                let mut iterations = Vec::with_capacity(pixel_count);
                let mut local_root_idx = Vec::with_capacity(pixel_count);
                let mut local_roots: Vec<Complex64> = Vec::new();

                for py in row_start..row_end {
                    for px in 0..WIDTH {
                        let x = X_MIN + (px as f64 / WIDTH as f64) * (X_MAX - X_MIN);
                        let y = Y_MAX - (py as f64 / HEIGHT as f64) * (Y_MAX - Y_MIN);
                        let z0 = Complex64::new(x, y);

                        match polynomial::iterate_newton(MAX_ITER, poly_ref, dpoly_ref, z0) {
                            Some((z_final, iter_count)) => {
                                let idx = match known_roots_ref {
                                    Some(known) => polynomial::nearest_known_root(known, z_final),
                                    None => polynomial::root_index(&mut local_roots, z_final),
                                };
                                iterations.push(iter_count);
                                local_root_idx.push(Some(idx));
                            }
                            None => {
                                iterations.push(MAX_ITER);
                                local_root_idx.push(None);
                            }
                        }
                    }
                }

                ThreadResult {
                    row_start,
                    iterations,
                    local_root_idx,
                    local_roots,
                }
            });
            handles.push(handle);
        }

        for handle in handles {
            raw_results.push(handle.join().unwrap());
        }
    });

    // Résolution des indices : selon si les racines étaient connues ou découvertes
    let resolution = if known_roots.is_some() {
        RootResolution::Known
    } else {
        let per_thread_roots: Vec<Vec<Complex64>> =
            raw_results.iter().map(|r| r.local_roots.clone()).collect();
        let (global_roots, mapping) = cluster_roots(&per_thread_roots, ROOT_MERGE_DIST);
        println!(
            "{} racines distinctes détectées (après clustering).",
            global_roots.len()
        );
        RootResolution::Discovered { mapping }
    };

    let mut img = RgbImage::new(WIDTH, HEIGHT);
    for (thread_idx, result) in raw_results.iter().enumerate() {
        for (i, (&iter_count, &local_idx)) in result
            .iterations
            .iter()
            .zip(result.local_root_idx.iter())
            .enumerate()
        {
            let px = (i as u32) % WIDTH;
            let py = result.row_start + (i as u32) / WIDTH;

            let color = match local_idx {
                Some(idx) => color_for(
                    resolution.global_index(thread_idx, idx),
                    iter_count,
                    poly.deg(),
                ),
                None => Rgb([0, 0, 0]),
            };
            img.put_pixel(px, py, color);
        }
    }

    img.save(&path).expect("échec de sauvegarde");
    println!("Image générée : {}", path);
}
