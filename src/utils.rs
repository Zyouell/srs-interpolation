//! This module contains utility functions we repeatedly use in the library.
use crate::error::InterpolationError;
use ark_ec::{
    short_weierstrass::{Affine, SWCurveConfig},
    CurveGroup,
};
use ark_ff::{Field, PrimeField};
use ark_std::{cfg_chunks_mut, vec::Vec, One, Zero};

use rayon::prelude::*;

/// This table is used as a lookup for swapping 8 bits into their reverse order.
const BYTE_SWAP_TABLE: [u8; 256] = [
    0, 128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240, 8, 136, 72, 200, 40,
    168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248, 4, 132, 68, 196, 36, 164, 100, 228, 20,
    148, 84, 212, 52, 180, 116, 244, 12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60,
    188, 124, 252, 2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242, 10, 138,
    74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250, 6, 134, 70, 198, 38, 166, 102,
    230, 22, 150, 86, 214, 54, 182, 118, 246, 14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94,
    222, 62, 190, 126, 254, 1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241,
    9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249, 5, 133, 69, 197, 37,
    165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245, 13, 141, 77, 205, 45, 173, 109, 237, 29,
    157, 93, 221, 61, 189, 125, 253, 3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179,
    115, 243, 11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251, 7, 135, 71,
    199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247, 15, 143, 79, 207, 47, 175, 111,
    239, 31, 159, 95, 223, 63, 191, 127, 255,
];

/// This function reverses the order of the `log_n` bits of `index` and returns the result.
/// The index is also assumed to be at most a [`u64`], although it will never be close to this large in practice.
pub(crate) fn bit_reverse(index: usize, log_n: usize) -> Result<usize, InterpolationError> {
    // Find out how many bytes are in our representation (rounded down)
    let num_bytes = log_n / 8;

    // We also need the number of bits in the last byte
    let leftover = log_n % 8;

    // Cast the usize to a slice of bytes
    let bytes = index.to_le_bytes();

    // Then we take the first `num_bytes` bytes and reverse them
    let mut reversed_bytes = Vec::<u8>::with_capacity(8);
    // We push the top (incomplete) byte first
    reversed_bytes.push(BYTE_SWAP_TABLE[bytes[num_bytes] as usize]);

    for i in 0..num_bytes {
        reversed_bytes.push(BYTE_SWAP_TABLE[bytes[num_bytes - 1 - i] as usize]);
    }
    reversed_bytes.resize(8, 0);

    let result = usize::from_le_bytes(reversed_bytes.as_slice().try_into().map_err(|_| {
        InterpolationError::InvalidParameters(format!(
            "Could not convert vector of length {} to array of length 8",
            reversed_bytes.len()
        ))
    })?);

    Ok(result >> (8 - leftover))
}

/// This function takes as input a mutable reference to a slice of affine points and a generator `g` as well as a round number.
/// It then mutates the slice in place to perform an FFT butterfly operation.
pub(crate) fn fft_round<E, F, const IS_FIRST_ROUND: bool>(
    points: &mut [Affine<E>],
    g: F,
    round_number: usize,
) -> Result<(), InterpolationError>
where
    E: SWCurveConfig<ScalarField = F>,
    F: PrimeField,
{
    let k = 1usize << round_number;
    let half = k >> 1;
    // If its the first round we don't have to rescale any points
    if !IS_FIRST_ROUND {
        cfg_chunks_mut!(points, k).try_for_each(|points: &mut [Affine<E>]| {
            distribute_powers(&mut points[half..], g);
            Result::<(), InterpolationError>::Ok(())
        })?;
    }
    let len = points.len() / k;

    let mut batch_inversion_accumulator = E::BaseField::one();
    let mut scratch_x: Vec<E::BaseField> = vec![E::BaseField::zero(); len * half];
    let mut scratch_y: Vec<E::BaseField> = vec![E::BaseField::zero(); len * half];
    points
        .chunks_mut(k)
        .enumerate()
        .for_each(|(i, points_chunk): (usize, &mut [Affine<E>])| {
            for j in 0..half {
                // We store the sum of the two x-coordinates in the scratch space
                scratch_x[half * i + j] += points_chunk[j].x + points_chunk[j + half].x;
                // Store y2 - y1 in the y scratch space
                scratch_y[half * i + j] += points_chunk[j + half].y - points_chunk[j].y;
                // Store x2 - x1 in the second points x-coordinate
                points_chunk[j + half].x -= points_chunk[j].x;
                // Store y2 + y1 in the second points y-coordinate
                points_chunk[j + half].y += points_chunk[j].y;
                // Multiply y2 + y1 by the product of the delta x's so far.
                points_chunk[j + half].y *= -batch_inversion_accumulator;
                // Multiply y2 - y1 by the product of the delta x's so far.
                scratch_y[half * i + j] *= batch_inversion_accumulator;
                // Update the accumulator with the denominator from this round.
                batch_inversion_accumulator *= points_chunk[j + half].x;
            }
        });

    batch_inversion_accumulator =
        batch_inversion_accumulator
            .inverse()
            .ok_or(InterpolationError::FieldError(
                "Could not invert batch inversion accumulator".to_string(),
            ))?;

    points.chunks_mut(k).enumerate().rev().for_each(
        |(i, points_chunk): (usize, &mut [Affine<E>])| {
            for j in (0..half).rev() {
                // Store (y2 + y1) / (x2 - x1) in the y-coordinate of the second point
                points_chunk[j + half].y *= batch_inversion_accumulator;
                // Store (y2 - y1) / (x2 - x1) in the y scratch space
                scratch_y[half * i + j] *= batch_inversion_accumulator;
                // Update the inversion accumulator
                batch_inversion_accumulator *= points_chunk[j + half].x;

                // Calculate new x and y coordinates for P - Q.
                points_chunk[j + half].x =
                    points_chunk[j + half].y * points_chunk[j + half].y - scratch_x[half * i + j];
                points_chunk[j + half].y = points_chunk[j + half].y
                    * (points_chunk[j].x - points_chunk[j + half].x)
                    - points_chunk[j].y;

                // Calculate new x and y coordinates for P + Q.
                let x = points_chunk[j].x;
                points_chunk[j].x =
                    scratch_y[half * i + j] * scratch_y[half * i + j] - scratch_x[half * i + j];
                points_chunk[j].y =
                    scratch_y[half * i + j] * (x - points_chunk[j].x) - points_chunk[j].y;
            }
        },
    );

    Ok(())
}

/// This function multiplies a list of affine points by the powers of `g` and stores the result in the same list.
pub(crate) fn distribute_powers<E, F>(coeffs: &mut [Affine<E>], g: F)
where
    E: SWCurveConfig<ScalarField = F>,
    F: PrimeField,
{
    let mut pow = g;
    coeffs.iter_mut().skip(1).for_each(|coeff| {
        *coeff = (*coeff * pow).into_affine();
        pow *= &g
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_reverse() {
        // Check that the bit reverse function is its own inverse
        let n = 1 << 10;
        for i in 0..n {
            let m = bit_reverse(i, 10).unwrap();
            assert_eq!(i, bit_reverse(m, 10).unwrap());
        }

        // Test against a known result
        let expected_result: [usize; 8] = [0, 4, 2, 6, 1, 5, 3, 7];
        for i in 0..8 {
            assert_eq!(expected_result[i], bit_reverse(i, 3).unwrap());
        }
    }
}
