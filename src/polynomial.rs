//! Ce module fournit les structures et fonctions pour représenter
//! et manipuler des polynômes à coefficients complexes.
use num_complex::Complex64;

#[derive(Debug)]
pub struct Polynome {
    coeffs: Vec<Complex64>,
}

impl Polynome {
    pub fn new(coeffs: Vec<Complex64>) -> Self {
        Polynome { coeffs }
    }

    pub fn from_roots(roots: &[Complex64]) -> Self {
        let mut coeffs = vec![Complex64::new(1.0, 0.0)];
        for &r in roots {
            // (z - r) s'écrit [ -r, 1 ] dans notre convention (indice = degré)
            coeffs = Self::multiply(&coeffs, &[-r, Complex64::new(1.0, 0.0)]);
        }
        Polynome { coeffs }
    }

    pub fn deg(&self) -> usize {
        self.coeffs.len() - 1
    }

    ///Evalue un polynome en un point z du plan complexe par la méthode d'Horner
    pub fn eval(&self, z: Complex64) -> Complex64 {
        self.coeffs
            .iter()
            .rev()
            .fold(Complex64::new(0.0, 0.0), |acc, c| acc * z + c)
    }

    ///Multiplie 2 polynomes entre eux
    fn multiply(a: &[Complex64], b: &[Complex64]) -> Vec<Complex64> {
        let mut result = vec![Complex64::new(0.0, 0.0); a.len() + b.len() - 1];
        for (i, &ca) in a.iter().enumerate() {
            for (j, &cb) in b.iter().enumerate() {
                result[i + j] += ca * cb;
            }
        }
        result
    }

    //Renvoie le polynome dérivé
    pub fn derivative(&self) -> Self {
        if self.deg() == 0 {
            return Polynome {
                coeffs: vec![Complex64::new(0., 0.)],
            };
        }

        let res = self.coeffs[1..]
            .iter()
            .enumerate()
            .map(|(i, &c)| c * Complex64::new((i + 1) as f64, 0.))
            .collect();

        Polynome { coeffs: res }
    }
}

//Applique l'algorithme de newton n fois sur z
pub fn iterate_newton(
    n: u32,
    p: &Polynome,
    dp: &Polynome,
    mut z: Complex64,
) -> Option<(Complex64, u32)> {
    for i in 0..n {
        let dz = dp.eval(z);

        if dz.norm() < 1e-12 {
            return Option::None;
        }
        let z_next = z - p.eval(z) / dz;

        if (z_next - z).norm() < 1e-9 {
            return Option::Some((z_next, i));
        }

        z = z_next;
    }

    None
}

pub fn root_index(found: &mut Vec<Complex64>, z: Complex64) -> usize {
    for (idx, r) in found.iter().enumerate() {
        if (z - r).norm() < 1e-4 {
            return idx;
        }
    }
    found.push(z);
    found.len() - 1
}

pub fn nearest_known_root(known: &[Complex64], z: Complex64) -> usize {
    known
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (z - **a).norm().partial_cmp(&(z - **b).norm()).unwrap())
        .map(|(idx, _)| idx)
        .unwrap()
}
