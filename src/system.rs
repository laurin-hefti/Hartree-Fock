use crate::orb;

pub struct Nucleus {
    pub pos: Vec<f64>,
    pub n: i32,
}

impl Nucleus {
    pub fn new(pos: Vec<f64>, n: i32) -> Nucleus{
        return Nucleus{pos:pos, n:n};
    }
}

pub struct Electron {
    pub pos: Vec<f64>,
    pub spin: f64,
    pub orb: Box<dyn Fn(&Vec<f64>, &Vec<f64>) -> f64+ 'static>,

    pub alpha: Vec<f64>,
    pub c_orb: Vec<f64>,

    pub num_basis: i32,
}

// alpha [l,m,n,alpha,N], first 3 elements are the angular momentum after the alpha value and then the norm

impl Electron {
    pub fn new(pos: Vec<f64>, spin: f64, 
        orb: Box<dyn Fn(&Vec<f64>, &Vec<f64>) -> f64+ 'static>, 
        alpha: Vec<f64>, c_orb: Vec<f64>) ->Electron {

            let num_basis: i32 = alpha.len() as i32 / 5;
            if num_basis != c_orb.len() as i32 {
                println!("error in electron : number of basis not correspont to number of coef!");
            }

        return Electron{pos:pos, spin:spin, orb:orb, alpha: alpha, c_orb: c_orb, num_basis: num_basis};
    }

    /*
    pub fn get_value(&self, c: &Vec<f64>) -> f64{
        let mut sum: f64= 0.0;
        for i in 0..self.num_basis {
            sum += self.c_orb[i as usize] * (self.orb)(c, &self.alpha[(i*5) as usize..(i*5 + 5) as usize].to_vec());
        }
        return sum; //normierung
    }
    */

    pub fn get_sing_value(&self, c: &Vec<f64>, i: i32) -> f64{
        return (self.orb)(c, &self.alpha[(i*5) as usize..(i*5 + 5) as usize].to_vec())
    }

    pub fn update_alpha(&mut self, a: f64, i: i32) {
        self.alpha[(i*5+3) as usize] = a;
        self.alpha[(i*5+4) as usize] = orb::NORM_GTO(self.alpha[(i*5+0) as usize], self.alpha[(i*5+1) as usize], self.alpha[(i*5+2) as usize], self.alpha[(i*5+3) as usize]);
    }

    pub fn print_alpha(&self){
        println!("{:?}", self.alpha);
    }
}

pub struct System {
    pub electrons : Vec<Electron>,
    pub nucleus : Vec<Nucleus>,

    pub num_basis: i32,
}

impl System {
    pub fn new(e1: Electron) -> System{
        let num_basis: i32 = e1.alpha.len() as i32 / 5;
        return System {
            electrons: vec![e1],
            nucleus: vec![],
            num_basis: num_basis,
        }
    }

    pub fn add_el(&mut self, e1: Electron){
        let num_basis: i32 = e1.alpha.len() as i32 / 5;
        self.electrons.push(e1);
        self.num_basis += num_basis;
    }

    pub fn add_nu(&mut self, nu: Nucleus){
        self.nucleus.push(nu);
    }
}
