use statrs::function::erf::erf;

const PI: f64= std::f64::consts::PI;
const e: f64= std::f64::consts::E;

// -- vector

pub fn len_vec(c: &Vec<f64>) -> f64{
    return (c[0].powf(2.0) + c[1].powf(2.0) + c[2].powf(2.0)).sqrt();
}

pub fn sub_vec(c1: &Vec<f64>, c2: &Vec<f64>) -> Vec<f64>{
    return vec![c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]];
}

pub fn add_vec(c1: &Vec<f64>, c2: &Vec<f64>) -> Vec<f64>{
    return vec![c1[0]+c2[0],c1[1]+c2[1],c1[2]+c2[2]];
}


pub fn scale_vec(c1: &Vec<f64>, s: f64) -> Vec<f64>{
    return vec![c1[0]*s,c1[1]*s,c1[2]*s];
}
// -- math

pub fn double_fac(n: i64) -> i64 {
    let mut res = 1;
    let mut k = n;

    while k > 1 {
        res *= k;
        k -=2;
    }

    return res;
}

pub fn boys(t: f64) -> f64{
    if t.abs() < 1e-8{
        return 1.0 - t/3.0 + t*t/10.0; // taylor expansion
    }

    if t > 1e6 {
        return 1.0 / (2.0*t);
    }

    return (1.0/2.0) * (PI/t).sqrt() * erf(t.sqrt() as f64) as f64;
}

// - linalg
#[derive(Clone)]
pub struct M {
    pub m: Vec<Vec<f64>>,
}

/*
[[0,0,0],
[0,0,0],
[0,0,0]]
*/

impl M {
    pub fn new(x: i32, y: i32)->M{
        return M{m:vec![vec![0.0; x as usize]; y as usize]}
    }

    pub fn print_m(&self){
        for v in &self.m{
            for el in v{
                print!("{} ", el);
            }
            println!("");
        }
    }

    pub fn add_m(m1: M, m2: M) -> M{
        if (m1.m.len() != m2.m.len() || m1.m[0].len() != m2.m[0].len()){
            println!("error in add Matrix: len of the matrix is not the same!");
        }

        let mut new_m = M::new(m1.m.len() as i32, m1.m[0].len() as i32);

        for i in 0..m1.m.len() {
            for y in 0..m1.m[0].len(){
                new_m.m[i][y] = m1.m[i][y] + m2.m[i][y];
            }
        }

        return new_m;
    }

    pub fn mult_m(m1: M, m2: M) -> M{
        if (m1.m.len() != m2.m.len() || m1.m[0].len() != m2.m[0].len()){
            println!("error in add Matrix: len of the matrix is not the same!");
        }

        let mut new_m = M::new(m1.m.len() as i32, m1.m[0].len() as i32);

        for i in 0..m1.m.len() {
            for y in 0..m1.m[0].len(){
                new_m.m[i][y] = m1.m[i][y] * m2.m[i][y];
            }
        }

        return new_m;
    }
}

