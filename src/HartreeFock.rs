use std::f32::INFINITY;

use crate::integrator;
use crate::system::*;
use crate::utils::*;
use nalgebra::DVector;
use nalgebra::{DMatrix, SymmetricEigen};

pub struct HF {
    pub S: M,
    pub Hcore: M,
    pub P: M,
    pub F: M,
    pub C: M,
    pub e: M,
}

impl HF {
    pub fn new(sys: &System) -> HF{
        let num_basis = sys.num_basis;
        return HF { S: M::new(num_basis, num_basis),
                    Hcore: M::new(num_basis, num_basis),
                    P: M::new(num_basis, num_basis),
                    F: M::new(num_basis, num_basis),
                    C: M::new(num_basis,num_basis),
                    e: M::new(num_basis, num_basis),
        }
    }

    pub fn calc_S(&mut self, sys: &System){
        for el in 0..sys.electrons.len() as i32{
            for el2 in 0..sys.electrons.len() as i32{
                for b1 in 0..sys.electrons[el as usize].num_basis{
                    for b2 in 0..sys.electrons[el2 as usize].num_basis{
                        
                        self.S.m[((el*sys.electrons[el as usize].num_basis)+b1) as usize][((el2*sys.electrons[el2 as usize].num_basis)+b2) as usize] = 
                        integrator::num_int_overlap(&sys.electrons[(el) as usize], 
                    &sys.electrons[(el2) as usize], b1, b2);
                    }
                }
            }
        }
    }

    pub fn calc_Hore(&mut self, sys: &System){
        for el in 0..sys.electrons.len() as i32{
            for el2 in 0..sys.electrons.len() as i32{
                for b1 in 0..sys.electrons[el as usize].num_basis{
                    for b2 in 0..sys.electrons[el2 as usize].num_basis{
                        
                        self.Hcore.m[((el*sys.electrons[el as usize].num_basis)+b1) as usize][((el2*sys.electrons[el2 as usize].num_basis)+b2) as usize] = 
                        integrator::num_int_pot_en(&sys.electrons[el as usize],
                             &sys.electrons[el2 as usize],  &sys.nucleus, b1, b2)
                             +
                             integrator::num_int_kin_en(&sys.electrons[el as usize], 
                                &sys.electrons[el2 as usize], b1, b2);
                    }
                }
            }
        }
    }

    
    pub fn calc_dens(&mut self, sys: &System){
        for v in 0..sys.num_basis {
            for y in 0..sys.num_basis{
                let mut sum =  0.0;
                for i in 0..sys.electrons.len()/2{
                    sum += self.C.m[v as usize][i] * self.C.m[y as usize][i];// * 2.0;
                }
                
                self.P.m[v as usize][y as usize] = 2.0 * sum;
            }
        }
    }

    pub fn clac_Fock(&mut self, sys: &System){
        for el in 0..sys.electrons.len() as i32{
            for el2 in 0..sys.electrons.len() as i32{
                for b1 in 0..sys.electrons[el as usize].num_basis{
                    for b2 in 0..sys.electrons[el2 as usize].num_basis{
                        
                        let mut sum_col = 0.0;
                        let mut sum_corr = 0.0;

                        for el3 in 0..sys.electrons.len() as i32{
                            for el4 in 0..sys.electrons.len() as i32{
                                for b3 in 0..sys.electrons[el3 as usize].num_basis{
                                    for b4 in 0..sys.electrons[el4 as usize].num_basis{
                                        
                                        let res:Vec<f32> = integrator::num_tow_el_int(&sys.electrons[el as usize],
                                            &sys.electrons[el2 as usize],
                                            &sys.electrons[el3 as usize],
                                            &sys.electrons[el4 as usize],
                                             b1,b2,b3,b4);

                                        sum_col += res[0] * self.P.m[((el3*sys.electrons[el3 as usize].num_basis)+b3) as usize][((el4*sys.electrons[el4 as usize].num_basis)+b4) as usize];
                                        sum_corr += res[1] * self.P.m[((el3*sys.electrons[el3 as usize].num_basis)+b3) as usize][((el4*sys.electrons[el4 as usize].num_basis)+b4) as usize];
                                    }
                                }
                            }
                        }
                        self.F.m[((el*sys.electrons[el as usize].num_basis)+b1) as usize][((el2*sys.electrons[el2 as usize].num_basis)+b2) as usize] = sum_col - 
                            0.5 * sum_corr + self.Hcore.m[((el*sys.electrons[el as usize].num_basis)+b1) as usize][((el2*sys.electrons[el2 as usize].num_basis)+b2) as usize];
                    }
                }
            }
        }
    }

    pub fn Roothaan_Haal(&mut self, sys: &System) -> (DMatrix<f32>, DVector<f32>){
        let o_flat: Vec<f32> = self.S.m.iter().flatten().cloned().collect();
        let overlap = DMatrix::from_row_slice(sys.num_basis as usize, sys.num_basis as usize, &o_flat[..]);
        let f_flat: Vec<f32> = self.F.m.iter().flatten().cloned().collect();
        let fock = DMatrix::from_row_slice(sys.num_basis as usize, sys.num_basis as usize, &f_flat[..]);

        let s_eigen = SymmetricEigen::new(overlap.clone());
        let s_vals = s_eigen.eigenvalues;
        let s_vecs = s_eigen.eigenvectors;

        let d_inv_sqrt = DMatrix::from_diagonal(
            &s_vals.map(|x|{
                if x > 0.0 { 1.0 / x.sqrt()}
                else {println!("error in diagonalisation? negativ value"); 0.0}
            })
        );
        let x = &s_vecs * d_inv_sqrt * s_vecs.transpose();

        let f_prime = x.transpose() * fock * &x;

        let f_eigen = SymmetricEigen::new(f_prime);
        let eps = f_eigen.eigenvalues;
        let c_tilde = f_eigen.eigenvectors;

        let c = &x * c_tilde;


        return (c, eps);
    }

    fn sort_MOs(c: &mut DMatrix<f32>, e: &mut DVector<f32>) {
        let mut e_copy = e.clone();
        let c_copy = c.clone();

        let mut min = e[0];
        let mut min_i = 0;
        for i in 0..e.len(){
            for j in 0..e.len(){
                if e_copy[j] < min {
                    min = e_copy[j];
                    min_i = j;
                }
            }

            e[i] = min;

            c.set_column(i, &c_copy.column(min_i as usize).clone());

            e_copy[min_i] = INFINITY;
            min = INFINITY;
        }
    }

    fn get_tot_E(&self, sys: &System) -> f32{
        let mut e_tot: f32 = 0.0;

        for v in 0..sys.num_basis {
            for y in 0..sys.num_basis {
                e_tot += self.P.m[v as usize][y as usize] * self.Hcore.m[v as usize][y as usize] + 0.5 * self.P.m[v as usize][y as usize] * self.F.m[v as usize][y as usize];
            }
        }

        return e_tot;
    }

    pub fn SCF(&mut self, sys: &System){
        println!("start scf");

        let DIIS: bool = true;

        self.calc_Hore(sys);
        self.calc_S(sys);
        self.calc_dens(sys);

        for it in 0..20 {
            println!("it : {}", it);
            self.clac_Fock(sys);

            let mut res: (DMatrix<f32>, DVector<f32>) = self.Roothaan_Haal(sys);

            HF::sort_MOs(&mut res.0, &mut res.1 );

            /*
            println!("P : ");
            self.P.print_m();
            println!("S : ");
            self.S.print_m();
            println!("F : ");
            self.F.print_m();
            */

            println!("F : ");
            self.F.print_m();

            let old_P: M = self.P.clone();

            self.calc_dens(sys);

            for v in 0..sys.num_basis{
                for y in 0..sys.num_basis{
                    self.C.m[v as usize][y as usize] = res.0[(v as usize,y as usize)];

                    self.P.m[v as usize][y as usize] *= 0.5;
                    self.P.m[v as usize][y as usize] += old_P.m[v as usize][y as usize]*0.5;
                }
            }
            
            /*
            println!(" en : \n{}", res.1);
            println!(" C :  \n{}", res.0);
            */

            println!(" C :  \n{}", res.0);

            println!("e_el : {}", self.get_tot_E(sys));
        }
    
        println!("end scf");
    }

}