pub const dx : f32 = 0.5;
pub const grid_dist : f32 = 2.0;    //volume to integrate
pub const int_b : i32 = (grid_dist / dx) as i32;

pub const dx_loss : f32 = 0.8;
pub const gird_dist_loss : f32 = 2.0;
pub const int_b_loss : i32 = (gird_dist_loss / dx_loss) as i32;

pub const nuc_radius_cut_off : f32 = 10e-10;
pub const orb_radius_cut_off : f32 = 10e-5;

pub const orb_val_cut_off : f32 = 0.01;