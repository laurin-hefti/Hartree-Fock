use crate::settings;
use crate::system::*;
use crate::utils::*;

pub const PI: f64= std::f64::consts::PI;
pub const e: f64= std::f64::consts::E;


pub fn num_int_pot_en(el: &Electron, el2: &Electron, nucleus: &Vec<Nucleus>, i: i32, j: i32) -> f64{
    let mut sum: f64= 0.0;

    let v_e1_e2: Vec<f64> = sub_vec(&el2.pos, &el.pos);

    for ix in -settings::int_b .. settings::int_b+1 {
        for iy in -settings::int_b .. settings::int_b+1 {
            for iz in -settings::int_b .. settings::int_b+1 {
                let x = ix as f64* settings::dx;
                let y = iy as f64* settings::dx;
                let z = iz as f64* settings::dx;

                let rel_cord = vec![x,y,z];

                for n in nucleus {
                    let r: f64= len_vec(&sub_vec(&add_vec(&n.pos, &rel_cord), &el.pos));

                    if r < settings::nuc_radius_cut_off{
                        continue;
                    }

                    sum += el.get_sing_value(&rel_cord, i) * 
                        el2.get_sing_value(&sub_vec(&rel_cord, &v_e1_e2), j)
                        * -n.n as f64/ r * settings::dx.powf(3.0);
                }
            }
        }
    }

    return sum;
}

pub fn num_int_kin_en(el: &Electron, el2: &Electron, i: i32, j: i32) -> f64{
    let mut sum: f64= 0.0;

    let v_e1_e2: Vec<f64> = sub_vec(&el2.pos, &el.pos);

    for ix in -settings::int_b .. settings::int_b+1 {
        for iy in -settings::int_b .. settings::int_b+1 {
            for iz in -settings::int_b .. settings::int_b+1 {
                let x = ix as f64* settings::dx;
                let y = iy as f64* settings::dx;
                let z = iz as f64* settings::dx;

                let x0 = vec![x,y,z];
                let x2 = vec![x+settings::dx,y,z];
                let x3 = vec![x-settings::dx,y,z];
                let y2 = vec![x,y+settings::dx,z];
                let y3 = vec![x,y-settings::dx,z];
                let z2 = vec![x,y,z+settings::dx];
                let z3 = vec![x,y,z-settings::dx];

                let v0 = el.get_sing_value(&x0, i);
                let vx1 = el.get_sing_value(&x2, i);
                let vx_1 = el.get_sing_value(&x3, i);
                let vy1 = el.get_sing_value(&y2, i);
                let vy_1 = el.get_sing_value(&y3, i);
                let vz1 = el.get_sing_value(&z2, i);
                let vz_1 = el.get_sing_value(&z3, i);

                let dx2 = (-2.0*v0 + vx1 + vx_1)/(settings::dx.powf(2.0));
                let dy2 = (-2.0*v0 + vy1 + vy_1)/(settings::dx.powf(2.0));
                let dz2 = (-2.0*v0 + vz1 + vz_1)/(settings::dx.powf(2.0));

                let laplace = dx2 + dy2 + dz2;

                sum += laplace * 
                el2.get_sing_value(&sub_vec(&vec![x,y,z], &v_e1_e2), j) * 
                settings::dx.powf(3.0);
            }
        }
    }

    return (-1.0/2.0) * sum;
}

pub fn num_tow_el_int(e1: &Electron, e2: &Electron, e3: &Electron, e4: &Electron,
    i:i32, j:i32, k:i32, l:i32) -> Vec<f64> {
    let mut sum: Vec<f64> = vec![0.0,0.0];

    let mut corr = false;
    if e1.spin == e2.spin {
        corr = true;
    }

    let v_e1_e2: Vec<f64> = sub_vec(&e2.pos, &e1.pos);

    for ix in -settings::int_b .. settings::int_b+1 {
        for iy in -settings::int_b .. settings::int_b+1 {
            for iz in -settings::int_b .. settings::int_b+1 {
                let x = ix as f64* settings::dx;
                let y = iy as f64* settings::dx;
                let z = iz as f64* settings::dx;

                let e1_cord = add_vec(&vec![x,y,z], &e1.pos);
                let e1_v = e1.get_sing_value(&vec![x,y,z], i);
                let e3_v = e3.get_sing_value(&vec![x,y,z], k);

                if e1_v.abs() < settings::orb_val_cut_off{
                    continue;
                }
                
                for i2x in -settings::int_b_loss .. settings::int_b_loss+1 {
                    for i2y in -settings::int_b_loss .. settings::int_b_loss+1 {
                        for i2z in -settings::int_b_loss .. settings::int_b_loss+1 {
                            let x2 = i2x as f64* settings::dx_loss;
                            let y2 = i2y as f64* settings::dx_loss;
                            let z2 = i2z as f64* settings::dx_loss;
                            
                            let e2_cord = add_vec(&vec![x2,y2,z2], &e2.pos);
                            let e2_v = e2.get_sing_value(&vec![x2,y2,z2], j);
                            let e4_v = e4.get_sing_value(&vec![x2,y2,z2], l);

                            let mut dist = len_vec(&sub_vec(&e1_cord, &e2_cord));

                            if dist < settings::orb_radius_cut_off{
                                dist = settings::orb_radius_cut_off;
                            }

                            sum[0] += e1_v*e3_v * (1.0/dist) * e2_v*e4_v * settings::dx.powf(3.0) * settings::dx_loss.powf(3.0);

                            if corr{
                                let e1_corr = e1.get_sing_value(&add_vec(&vec![x,y,z], &v_e1_e2), i);
                                let e2_corr = e2.get_sing_value(&sub_vec(&vec![x2,y2,z2,], &v_e1_e2), j);

                                sum[1] += (e3_v * e1_corr) * (1.0/dist) * (e4_v * e2_corr) * settings::dx.powf(3.0) * settings::dx_loss.powf(3.0);
                            }
                        }
                    }
                }      
            }
        }
    }

    //sum[0] *= 1.0;
    //sum[1] *= -0.5;
    return sum;
}

pub fn num_int_overlap(e1: &Electron, e2: &Electron, i:i32, j:i32) -> f64{
    let mut sum: f64= 0.0;

    let v_e1_e2: Vec<f64> = sub_vec(&e2.pos, &e1.pos);

    for ix in -settings::int_b .. settings::int_b+1 {
        for iy in -settings::int_b .. settings::int_b+1 {
            for iz in -settings::int_b .. settings::int_b+1 {
                let x = ix as f64* settings::dx;
                let y = iy as f64* settings::dx;
                let z = iz as f64* settings::dx;

                sum += e1.get_sing_value(&vec![x,y,z], i) * 
                e2.get_sing_value(&sub_vec(&vec![x,y,z], &v_e1_e2),j) *
                settings::dx.powf(3.0);
            }
        }
    }

    return sum;
}

//for primitives
pub fn an_s_overlap(e1: &Electron, e2: &Electron, i: i32, j: i32) -> f64{
    let a: f64= e1.alpha[(i*5+3) as usize];
    let b: f64= e2.alpha[(j*5+3) as usize];
    let n1: f64 = e1.alpha[(i*5+4) as usize];
    let n2: f64 = e1.alpha[(j*5+4) as usize];

    return n1 * n2 * (PI / (a+b)).powf(3.0/2.0) * 
        e.powf(-((a*b)/(a+b)) * len_vec(&sub_vec(&e1.pos, &e2.pos)).powf(2.0));
}

pub fn an_s_kin(e1: &Electron, e2: &Electron, i: i32, j: i32) -> f64{
    let a: f64= e1.alpha[(i*5+3) as usize];
    let b: f64= e2.alpha[(j*5+3) as usize];
    let n1: f64 = e1.alpha[(i*5+4) as usize];
    let n2: f64 = e1.alpha[(j*5+4) as usize];
    let y: f64= (a*b)/(a+b);

    let s: f64= n1 * n2 * (PI / (a+b)).powf(3.0/2.0) * 
        e.powf(-y * len_vec(&sub_vec(&e1.pos, &e2.pos)).powf(2.0));

    return (y * (3.0-2.0*y*len_vec(&sub_vec(&e1.pos, &e2.pos)).powf(2.0)) * s);
}

pub fn an_s_pot(e1: &Electron, e2: &Electron, nucleus: &Vec<Nucleus>, i: i32, j: i32) -> f64{
    let a: f64= e1.alpha[(i*5+3) as usize];
    let b: f64= e2.alpha[(j*5+3) as usize];
    let n1: f64 = e1.alpha[(i*5+4) as usize];
    let n2: f64 = e1.alpha[(j*5+4) as usize];
    let y: f64= (a*b)/(a+b);

    let p_center: Vec<f64> = scale_vec(
        &add_vec(&scale_vec(&e1.pos, a), &scale_vec(&e2.pos, b)),
    1.0/(a+b));

    let mut sum: f64= 0.0;

    for z in 0..nucleus.len(){
        sum += - nucleus[z].n as f64 * n1 * n2 *((2.0*PI)/(a+b)) *
            e.powf(-y * len_vec(&sub_vec(&e1.pos, &e2.pos)).powf(2.0)) *
            boys((a*b) * len_vec(&sub_vec(&p_center, &nucleus[z].pos)));
    }

    return sum;
}

//for primitives
pub fn an_s_two_el(e1: &Electron, e2: &Electron, e3: &Electron, e4: &Electron,
    i:i32, j:i32, k:i32, l:i32) -> Vec<f64> {
    let mut sum: Vec<f64> = vec![0.0,0.0];

    let mut corr = false;

    if e1.spin == e2.spin {
        corr = true;
    }

    let  a: f64= e1.alpha[(i*5+3) as usize];
    let b: f64= e2.alpha[(j*5+3) as usize];
    let y: f64= e3.alpha[(k*5+3) as usize];
    let d: f64= e4.alpha[(l*5+3) as usize];
    let n1: f64 = e1.alpha[(i*5+4) as usize];
    let n2: f64 = e1.alpha[(j*5+4) as usize];
    let n3: f64 = e1.alpha[(k*5+4) as usize];
    let n4: f64 = e1.alpha[(l*5+4) as usize];
    let mut p: f64= a+b;
    let mut q: f64= y+d;
    if p < 10e-6 {p = 10e-6;}       //not necesary because alphy is not so little
    if q < 10e-6 {q = 10e-6;}
    let mut pq = p+q;
    let yab: f64= (a*b)/p;
    let yyd: f64= (y*d)/q;

    let p_center: Vec<f64> = scale_vec(
        &add_vec(&scale_vec(&e1.pos, a), &scale_vec(&e2.pos, b)),
        1.0/p);

    let q_center: Vec<f64> = scale_vec(
        &add_vec(&scale_vec(&e3.pos, y), &scale_vec(&e4.pos, d)),
        1.0/q);

    let mut v_e1_e2: f64 = len_vec(&sub_vec(&e1.pos, &e2.pos)).powf(2.0);
    let mut v_e3_e4: f64 = len_vec(&sub_vec(&e3.pos, &e4.pos)).powf(2.0);
    let mut v_c1_c2: f64 = len_vec(&sub_vec(&p_center, &q_center)).powf(2.0);

    if v_e1_e2 < 10e-6 {v_e1_e2 = 10e-6;}       //only check one side because is power of 2
    if v_e3_e4 < 10e-6 {v_e3_e4 = 10e-6;}
    if v_c1_c2 < 10e-6 {v_c1_c2 = 10e-6;}

    sum[0] = n1 * n2 * n3 * n4 *  (((2.0*PI).powf(5.0/2.0)) / ((p*q)*(pq).sqrt())) *
        e.powf(-yab * v_e1_e2
                - yyd * v_e3_e4)
            * boys((p*q)/(pq) * v_c1_c2);

    if corr {
        let mut p_cor = a+d;
        let mut q_cor = y+b;
        if p_cor < 10e-6 {p_cor = 10e-6;}
        if q_cor < 10e-6 {q_cor = 10e-6;}
        let pq_cor = p_cor + q_cor;
        let yad_cor = (a*d) / p_cor;
        let yyb_cor = (y*b) / q_cor;
       
        let p_center_cor: Vec<f64> = scale_vec(
        &add_vec(&scale_vec(&e1.pos, a), &scale_vec(&e4.pos, d)),
        1.0/p_cor);

        let q_center_cor: Vec<f64> = scale_vec(
        &add_vec(&scale_vec(&e3.pos, y), &scale_vec(&e2.pos, b)),
        1.0/q_cor);

        let mut v_e1_e4: f64 = len_vec(&sub_vec(&e1.pos, &e4.pos)).powf(2.0);
        let mut v_e2_e3: f64 = len_vec(&sub_vec(&e3.pos, &e2.pos)).powf(2.0);
        let mut v_c1_c2_corr: f64 = len_vec(&sub_vec(&p_center_cor, &q_center_cor)).powf(2.0);
        
        if v_e1_e4 < 10e-6 {v_e1_e4 = 10e-6;}
        if v_e2_e3 < 10e-6 {v_e2_e3 = 10e-6;}
        if v_c1_c2_corr < 10e-6 {v_c1_c2_corr = 10e-6;}
 
        sum[1] = n1 * n2 * n3 * n4 * (((2.0*PI).powf(5.0/2.0)) / ((p_cor*q_cor)*(pq_cor).sqrt())) *
        e.powf(-yad_cor * v_e1_e4
                - yyb_cor * v_e2_e3)
            * boys((p_cor*q_cor)/(pq_cor) * v_c1_c2_corr);
    }
    
    return sum;
    }