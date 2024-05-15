use ark_ec::{
    short_weierstrass::{Affine, SWCurveConfig},
    CurveGroup,
};

use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{cfg_into_iter, cfg_iter};
use error::InterpolationError;
use rayon::prelude::*;
use utils::{bit_reverse, fft_round};
mod error;
mod utils;
// To begin with we assume our SRS is given to us in ascending order of the powers of tau,
// that is [1], [tau], [tau^2], [tau^3], ... [tau^d]. Since we are doing many point additions
// we want to use the batch inversion trick and work in affine form to speed this up.

/// This function performs batch affine addition on a series of pairs of affine points in short Weierstrass form.
/// We assume the points have been arranged in a manner so that points to be added are adjacent to one another.
pub fn srs_to_lagrange<E, F>(points: &[Affine<E>]) -> Result<Vec<Affine<E>>, InterpolationError>
where
    E: SWCurveConfig<ScalarField = F>,
    F: PrimeField,
{
    // First we check that the number of points is a power of two.
    let log_point_size = points.len().ilog2() as usize;
    let point_size = 1usize << log_point_size;
    if points.len() != point_size {
        return Err(InterpolationError::SizeError);
    }

    if points.len() == 1 {
        return Ok(points.to_vec());
    }

    // First we order the points so that it is convenient to perform the FFT style operation.
    let mut ordered_points = cfg_into_iter!(0..points.len())
        .map(|i| Ok(points[bit_reverse(i, log_point_size)?]))
        .collect::<Result<Vec<Affine<E>>, InterpolationError>>()?;

    // We need the |2^log_point_size|th root of unity in the field.
    let domain =
        Radix2EvaluationDomain::<F>::new(point_size).ok_or(InterpolationError::SizeError)?;
    let gen = domain.group_gen_inv();

    // Then we perform the FFT style operation.
    for i in 1..=log_point_size {
        // In each round we take the point_size >> i th root of unity
        let prim_root = gen.pow(&[(point_size >> i) as u64]);
        if i != 1 {
            fft_round::<E, F, false>(&mut ordered_points, prim_root, i)?;
        } else {
            fft_round::<E, F, true>(&mut ordered_points, prim_root, i)?;
        }
    }

    // Finally we rescale all the points by the size of the domain
    let domain_size_inv =
        domain
            .size_as_field_element
            .inverse()
            .ok_or(InterpolationError::FieldError(
                "Could not invert domain size".to_string(),
            ))?;
    Ok(cfg_iter!(ordered_points)
        .map(|point| (*point * domain_size_inv).into_affine())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{g1::Config as BnConfig, Fr};
    use ark_ec::{
        scalar_mul::fixed_base::FixedBase, short_weierstrass::Projective, CurveGroup,
        VariableBaseMSM,
    };

    use ark_std::{rand::RngCore, One, UniformRand};

    /// We only generate the G1 points for testing purposes.
    pub(crate) fn gen_srs_for_testing<E: SWCurveConfig, R: RngCore>(
        rng: &mut R,
        max_degree: usize,
    ) -> Result<Vec<Affine<E>>, InterpolationError> {
        let beta = E::ScalarField::rand(rng);
        let g = Projective::<E>::rand(rng);

        let mut powers_of_beta = vec![E::ScalarField::one()];

        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBase::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g =
            FixedBase::msm::<Projective<E>>(scalar_bits, window_size, &g_table, &powers_of_beta);

        let powers_of_g = Projective::<E>::normalize_batch(&powers_of_g);

        Ok(powers_of_g)
    }

    #[test]
    fn test_srs_interpolation() -> Result<(), InterpolationError> {
        test_srs_interpolation_helper::<BnConfig, Fr>();
        Ok(())
    }

    fn test_srs_interpolation_helper<E, F>()
    where
        E: SWCurveConfig<ScalarField = F>,
        F: PrimeField,
    {
        let rng = &mut ark_std::test_rng();
        for i in 5..10 {
            let max_degree = (1usize << i) - 1;

            let domain = Radix2EvaluationDomain::<F>::new(max_degree + 1).unwrap();
            let srs = gen_srs_for_testing::<E, _>(rng, max_degree).unwrap();

            let evals = (0..(max_degree + 1))
                .map(|_| F::rand(rng))
                .collect::<Vec<F>>();

            let coeffs = domain.ifft(&evals);

            let coeff_commitment = Projective::<E>::msm_bigint(
                &srs,
                &coeffs.iter().map(|x| x.into_bigint()).collect::<Vec<_>>(),
            )
            .into_affine();

            let lagrange_srs = srs_to_lagrange::<E, F>(&srs).unwrap();

            let lagrange_commitment = Projective::<E>::msm_bigint(
                &lagrange_srs,
                &evals.iter().map(|x| x.into_bigint()).collect::<Vec<_>>(),
            )
            .into_affine();

            assert_eq!(coeff_commitment, lagrange_commitment);
        }
    }
}
