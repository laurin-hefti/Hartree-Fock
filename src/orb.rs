use crate::utils::*;

pub const PI: f32 = std::f32::consts::PI;
pub const e: f32 = std::f32::consts::E;

// alpha [l,m,n,alpha,N], first 3 elements are the angular momentum after the alpha value and then the norm

pub fn H1S_rel(c: &Vec<f32>, alpha: &Vec<f32>) -> f32  {
    return (1.0/PI.sqrt())*e.powf(-len_vec(c));
}

pub fn GTO(c: &Vec<f32>, alpha: &Vec<f32>) -> f32{
    return alpha[4] * c[0].powf(alpha[0]) * c[1].powf(alpha[1]) * c[2].powf(alpha[2]) * e.powf(-alpha[3] * len_vec(c).powf(2.0));
}

pub fn NORM_1S_GTO(a: f32) -> f32{
    return (2.0*a / PI).powf(3.0/4.0);
}

pub fn NORM_GTO(l: f32, m: f32, n: f32, a: f32) -> f32{
    return (
        (
            2.0_f32.powf(l+m+n+(3.0/2.0)) * a.powf(l+m+n+3.0/2.0)
        ) /
        (
            PI.powf(3.0/2.0) * double_fac((2.0*l-1.0)as i64) as f32 *
            double_fac((2.0*m-1.0)as i64) as f32 *
            double_fac((2.0*n-1.0)as i64) as f32
        )
    ) .powf(1.0/2.0)
}

pub fn STO(){}
