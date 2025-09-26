pub mod system;
pub mod orb;
pub mod integrator;
pub mod settings;
pub mod utils;
pub mod HartreeFock;

use crate::system::*;
use crate::HartreeFock::*;

fn main() {
    let el: system::Electron = Electron::new(vec![0.0,0.0,0.0], 
        -1.0/2.0, 
        Box::new(orb::GTO),
        vec![

        0.0,0.0,0.0,1.2,orb::NORM_GTO(0.0, 0.0, 0.0, 1.2),
        0.0,0.0,0.0,0.8,orb::NORM_GTO(0.0, 0.0, 0.0, 0.8),
        
        0.0,0.0,0.0,0.3,orb::NORM_GTO(0.0, 0.0, 0.0, 0.3),
        0.0,0.0,0.0,0.2,orb::NORM_GTO(0.0,0.0,0.0,0.2)],
        vec![0.2,0.2,0.2,0.2,0.2, 0.2],
    );
    let nu: system::Nucleus = Nucleus{n: 1, pos:vec![0.0,0.0,0.0]};

    let el2: system::Electron = Electron::new(vec![1.4,0.0,0.0], 
        1.0/2.0, 
        Box::new(orb::GTO),
        vec![

        0.0,0.0,0.0,1.2,orb::NORM_GTO(0.0, 0.0, 0.0, 1.2),
        0.0,0.0,0.0,0.8,orb::NORM_GTO(0.0, 0.0, 0.0, 0.8), 

        0.0,0.0,0.0,0.3,orb::NORM_GTO(0.0, 0.0, 0.0, 0.3),
        0.0,0.0,0.0,0.2,orb::NORM_GTO(0.0, 0.0, 0.0, 0.2)],
        vec![0.2,0.2,0.2,0.2, 0.2, 0.2],
    );
    let nu2: system::Nucleus = Nucleus{n: 1, pos:vec![1.4,0.0,0.0]};
    
    let mut sys: System = System::new(el);
    sys.add_el(el2);
    sys.add_nu(nu);
    sys.add_nu(nu2);
    

    let mut hf1: HF = HF::new(&sys);
    hf1.SCF(&sys);
}
