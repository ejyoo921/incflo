#include <incflo.H>

using namespace amrex;

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

                //TODO: Here, we need to update cp et al., based on temperature
                // Is our tracer always temperature? 
                cp_fe = 300;
                cp_slg = 1000;
                cp_arr(i,j,k,0) = cp_fe; //TODO: define ID = 0 is cp_fe
                cp_arr(i,j,k,1) = cp_slg; //TODO: define ID = 1 is cp_slg
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

// void incflo::update_properties(int lev, MultiFab& t_prop)
// {

//     // for (MFIter mfi(t_prop); mfi.isValid(); ++mfi)
//     // {
//     //     Array4<Real> t_prop_arr = t_prop.array(mfi);

//     //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//     //         incflo::update_thermal_properties_and_phases(i, j, k, t_prop_arr);
//     //     });
//     // }
// }