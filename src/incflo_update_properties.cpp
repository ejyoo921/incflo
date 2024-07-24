#include <incflo.H>

using namespace amrex;

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real bound01(const Real q0)
{   //bounds a quantity between 0 and 1
    amrex::Real q1 = q0;
    if (q0<0.0)
    {
        q1 = 0.0;
    }
    if (q0 > 1.0)
    {
        q1 = 1.0;
    }
        return q1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cp_1(const amrex::Real Temp)
{   //cp calc for solid pellet phase 1  
    amrex::Real cp_1a =0.4707;    
    amrex::Real T1 =  600.0;
    amrex::Real T2 = 1027.5;
    amrex::Real T2A= T1;
    amrex::Real T3 = 1227.5;
    amrex::Real T3A= T2;
    amrex::Real T4 =  1809.5;
    amrex::Real T4A=  T3;
    amrex::Real T5 = 1818.0;
    amrex::Real T5A =T4;
    amrex::Real cp1 = cp_1a;
    amrex::Real cp2 = cp1;
    amrex::Real cp3 = 1.0078;
    amrex::Real cp4 = 0.6066;
    amrex::Real cp5 = 0.7442;
    amrex::Real cp6 = 0.8951;
    amrex::Real a2 = 0.001256;
    amrex::Real a3 =-0.0020085;
    amrex::Real a4 = 0.0002362;
    amrex::Real a5 = 0.1352;
    amrex::Real c3 = 0.78478;
    amrex::Real c4 = 0.43746;
    amrex::Real c5 =11.331;
    amrex::Real o3 = 7.4375;
    amrex::Real o4 = 11.6875;
    amrex::Real o5 = 8.5;
    amrex::Real Tm3 = 1187.9375;
    amrex::Real Tm4 = 1670.35;
    amrex::Real Tm5= T4;

    if (Temp<= T1)
    { cp_1a = cp1;
    }
    else
    { if (Temp <= T2)
        { cp_1a = cp2 + a2*(Temp-T2A);
        }
        else
        { if (Temp <= T3)
        { cp_1a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o3),2));
            }
            else
        { if (Temp <= T4)
            {cp_1a = cp4 + a4*(Temp-T4A) + c4 * exp(-0.5*pow(((Temp-Tm4)/o4),2)) + c5 * exp(-0.5*pow(((Temp-Tm5)/o5),2));
            }
        else
        { if (Temp <= T5)
                {cp_1a = cp5 + a5*(Temp-T5A) + c5 * exp(-0.5*pow(((Temp-Tm5)/o5),2));
                }
            else
        {cp_1a = cp6+ c5 * exp(-0.5*pow(((Temp-Tm5)/o5),2));
        }}}}}
    
    return 1000.0*cp_1a;
}

AMREX_GPU_DEVICE AMREX_INLINE
amrex::Real cp_2(const amrex::Real Temp)
{   //cp calc for solid pellet phase 2  
    amrex::Real cp_2a = 0.5867;
    
    amrex::Real T1 = 678.875;
    amrex::Real T2 = 997.570;
    amrex::Real T2A= T1;
    amrex::Real T3 = 1040.25;
    amrex::Real T3A= 1016.875;
    amrex::Real T4 = 2073.0;
    amrex::Real T4A= 1059.4;
    amrex::Real cp1 = 0.5867;
    amrex::Real cp2 = cp1;
    amrex::Real cp3 = 1.1242;
    amrex::Real cp4 = 0.8144;
    amrex::Real cp5 = 0.9650;
    amrex::Real a2 = 0.000167513;
    amrex::Real a3 =-0.00729;
    amrex::Real a4 = 0.00014857;
    amrex::Real c2 = 2.5489;
    amrex::Real c3 = c2;
    amrex::Real c4 = 8.6987;
    amrex::Real o2 = 9.5625;
    amrex::Real o3 = o2;
    amrex::Real o4 = 10.625;
    amrex::Real Tm2= 1007.3125;
    amrex::Real Tm3 = Tm2;
    amrex::Real Tm4 = 1429.125;

    if (Temp<= T1)
{ cp_2a = cp1;
}
    else
{ if (Temp <= T2)
    { cp_2a = cp2 + a2*(Temp-T2A) + c2 * exp(-0.5*pow(((Temp-Tm2)/o2),2));
    }
        else
    { if (Temp <= T3)
        { cp_2a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o3),2));
        }
            else
        { if (Temp <= T4)
            {  cp_2a = cp4 + a4*(Temp-T4A) + c4 * exp(-0.5*pow(((Temp-Tm4)/o4),2));
            }
                else
            {cp_2a = cp5;
            }}}}

    return 1000.0*cp_2a;
}

AMREX_GPU_DEVICE AMREX_INLINE
amrex::Real cp_5(const amrex::Real Temp)
{ //cp calc for slag
    amrex::Real cp_5a = 1.0824;

    amrex::Real T1 =  834.125;
    amrex::Real T2 = 1053.0 ;
    amrex::Real T2A= T1;
    amrex::Real T3 = 1478.0;
    amrex::Real T3A= T2;
    amrex::Real T4 = 2073.0;
    amrex::Real T4A= T3;
    amrex::Real cp1 = 1.0824;
    amrex::Real cp2 = cp1;
    amrex::Real cp3 = 1.1844;
    amrex::Real cp4 = cp3;
    amrex::Real cp5 = 1.1794;
    amrex::Real a2 = 0.0004277;
    amrex::Real a3 =-0.0001121;
    amrex::Real a4 = 0.0006715;
    amrex::Real c1 = 0.4370;
    amrex::Real c2 = c1;
    amrex::Real c3 = 4.4169;
    amrex::Real c4 = c3;
    amrex::Real o1 = 59.5;
    amrex::Real o2 = o1;
    amrex::Real o3 = 48.3475;
    amrex::Real o4 = o3;
    amrex::Real Tm1= 893.625;
    amrex::Real Tm2= Tm1;
    amrex::Real Tm3 = 1574.6875;
    amrex::Real Tm4 = Tm3;

    if (Temp<= T1)
    { cp_5a = cp1 + c1 * exp(-0.5*pow(((Temp-Tm1)/o1),2));;
    }
    else
    { if (Temp <= T2)
    { cp_5a = cp2 + a2*(Temp-T2A) + c2 * exp(-0.5*pow(((Temp-Tm2)/o2),2));
    }
        else
    { if (Temp <= T3)
        { cp_5a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o3),2));
        }
            else
        { if (Temp <= T4)
            {  cp_5a = cp4 + a4*(Temp-T4A) + c4 * exp(-0.5*pow(((Temp-Tm4)/o4),2));
            }
                else
            {cp_5a = cp5;
            }}}}

    return 1000.0*cp_5a;
}

amrex::Real compute_cp(const amrex::Real Temp, int phase)
{
    // calculate cp based on phase and temp 
    // unit: J/K/kg
    Real cpcalc = 1000.0;
    switch (phase) {
        case 1:
           cpcalc = cp_1(Temp);
        break;
        case 2:
           cpcalc = cp_2(Temp);
        break;
        case 3:
           cpcalc = 2000.0; //cp calc for 1 % C bath
        break;
        case 4:
           cpcalc = 2000.0; //cp calc for 3% carbon bath  
        break;
        case 5:
           cpcalc = cp_5(Temp);
        break;
      }
    return cpcalc;
} 
}


void incflo::update_properties ()
{
    BL_PROFILE("incflo::update_properties");

    // update thermal propeties
    for(int lev = 0; lev <= finest_level; ++lev) {
        // you need to deal with ghost cells
        int ng = nghost_state();
        fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng); 

        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();
        
        amrex::ParmParse pp("prob");
        // Extract position and velocities
        amrex::Vector<amrex::Real> rads;
        amrex::Vector<amrex::Real> centx;
        amrex::Vector<amrex::Real> centy;
        amrex::Vector<amrex::Real> centz;
        int npellets = 0;

        pp.getarr("pellet_rads",  rads);    
        pp.getarr("pellet_centx", centx);
        pp.getarr("pellet_centy", centy);
        pp.getarr("pellet_centz", centz);
        pp.get("npellets", npellets);

        auto& ld = *m_leveldata[lev];
        for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Box const& bx = mfi.tilebox();
            Box const& bx = mfi.growntilebox(); //bigger box (with ghosts)
            Array4<Real> temp_arr = ld.tracer.array(mfi);
            Array4<Real> cp_arr = ld.cp.array(mfi); 

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // incflo::update_thermal_properties_and_phases(i, j, k, t_prop_arr, prob_lo, prob_hi, dx);
                // For SteelSlag case 
                amrex::Real Temp = temp_arr(i, j, k); // Need to bring actual Temp vals
                amrex::Real cp_fe,  cond_fe,  dens_fe; // Fe
                amrex::Real cp_slg, cond_slg, dens_slg; // Slag
                amrex::Real sol_fe, mol_fe, sol_slg, mol_slg; // Phases
                amrex::Real vfrac_fe;

                //TODO: Here, we need to update cp et al., based on temperature
                // Is our tracer always temperature? 
                cp_fe   = compute_cp(Temp, 2);
                cp_slg  = compute_cp(Temp, 5);
                // Q: how can I determine where to put each cp?
                vfrac_fe = bound01((temp_arr(i,j,k)-dens_slg)/(dens_fe-dens_slg));
                cp_arr(i,j,k) = cp_slg*(1.0-vfrac_fe)+cp_fe*vfrac_fe;
            });
        }
    }

    // set zero velocity - inside pellet
    ParmParse pp("incflo");
    bool m_zero_vel = false;
    pp.queryAdd("zero_vel", m_zero_vel);

    if (m_zero_vel) // Make zero velocity 
    {
        amrex::Print() << "Let's impose zero velocity" << "\n";
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            auto const& dx = geom[lev].CellSizeArray();
            amrex::ParmParse pp("prob");
            // Extract position and velocities
            amrex::Vector<amrex::Real> rads;
            amrex::Vector<amrex::Real> centx;
            amrex::Vector<amrex::Real> centy;
            amrex::Vector<amrex::Real> centz;

            int npellets = 0;
            Real density_p = 1.;

            pp.get("npellets", npellets);
            pp.get("density_p", density_p);

            pp.getarr("pellet_rads",  rads);    
            pp.getarr("pellet_centx", centx);
            pp.getarr("pellet_centy", centy);
            pp.getarr("pellet_centz", centz);

            auto& ld = *m_leveldata[lev];
            for (int lev = 0; lev <= finest_level; lev++)
            {
                auto const& problo = geom[lev].ProbLoArray();
                for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Array4<Real> const& vel = ld.velocity.array(mfi);
                    Box const& bx = mfi.tilebox();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real x = problo[0] + Real(i+0.5)*dx[0];
                        Real y = problo[1] + Real(j+0.5)*dx[1];
                        Real z = problo[2] + Real(k+0.5)*dx[2];

                        int inside_pellet = 0;
                        for(int np = 0; np < npellets; np++)
                        {

                            Real dist2 = std::pow(x - centx[np], 2.0)+                
                                        std::pow(y - centy[np], 2.0)+                
                                        std::pow(z - centz[np], 2.0); 
                            
                            if(dist2 < std::pow(rads[np], 2.0))
                            {              
                                inside_pellet = 1;
                                break;
                            }
                        }
                        if (inside_pellet)
                        {   // no internal flow
                            vel(i,j,k,0) = Real(0.0);
                            vel(i,j,k,1) = Real(0.0);
                            vel(i,j,k,2) = Real(0.0);
                        }
                    });
                }
            }
        }
    }
}

