#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    int npellets=1;
    amrex::Real* pellet_rads = nullptr;
    amrex::Real* pellet_centx = nullptr;
    amrex::Real* pellet_centy = nullptr;
    amrex::Real* pellet_centz = nullptr;
                                 //
    amrex::Real rho_liq_fe = 7500.0; //kg/m3
    amrex::Real cp_liq_fe  = 680.0;  //J/K/kg
    amrex::Real k_liq_fe   = 40.0;  //W/m/K
    amrex::Real rho_solid_fe = 7500.0;  //kg/m3
    amrex::Real cp_solid_fe  = 745.0;   //J/K/kg
    amrex::Real k_solid_fe   = 4.5;     //W/m/K
    amrex::Real Tmelt_fe = 1450.0+273.0; //K
    amrex::Real Tmelt1_fe=1460.0+273.0;
    amrex::Real Tmelt2_fe=1480.0+273.0;
    amrex::Real lheat_fe=270000.0;
    amrex::Real Cpmelt_fe=cp_liq_fe+lheat_fe/(Tmelt2_fe-Tmelt1_fe);
    
    amrex::Real rho_liq_slg = 2550.0; //kg/m3
    amrex::Real cp_liq_slg  = 1004.0;  //J/K/kg
    amrex::Real k_liq_slg   = 1.1715;  //W/m/K
    amrex::Real rho_solid_slg = 2550.0;  //kg/m3
    amrex::Real cp_solid_slg  = 1004.0;   //J/K/kg
    amrex::Real k_solid_slg   =  1.0;     //W/m/K
    amrex::Real Tmelt_slg = 1594.0; //K
    amrex::Real Tmelt1_slg=1589.0;
    amrex::Real Tmelt2_slg=1609.0;
    amrex::Real lheat_slg=9205904.0;
    amrex::Real Cpmelt_slg=cp_liq_slg+lheat_slg/(Tmelt2_slg-Tmelt1_slg);
    
    amrex::Real Tinit_solid_fe = 80.0+273.0; //K
    amrex::Real Tinit_liq_slg = 1540.0+273.0; //K

};
#endif
