#include <incflo.H>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = static_cast<Real>(ParallelDescriptor::second());

    // Compute time step size
    int initialisation = ( m_dt < 0 );
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    // Set new and old time to correctly use in fillpatching
    for(int lev = 0; lev <= finest_level; lev++)
    {
        amrex::MultiFab::Copy(t_prop_o[lev], t_prop[lev], 
                                  0, 0, t_prop[lev].nComp(), 0);
        m_t_old[lev] = m_cur_time;
        m_t_new[lev] = m_cur_time + m_dt;
    }

    if (m_verbose > 0)
    {
        amrex::Print() << "\nStep " << m_nstep + 1
                       << ": from old_time " << m_cur_time
                       << " to new time " << m_cur_time + m_dt
                       << " with dt = " << m_dt << ".\n" << std::endl;
    }

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();
    //EY
    // copy_from_new_to_old_t_prop();

    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }

#ifdef AMREX_USE_EB
    if (m_eb_flow.enabled) {
        amrex::Print() << "EB flow enabled" << "\n";  //EY
       for (int lev = 0; lev <= finest_level; ++lev) {
         set_eb_velocity(lev, m_t_old[lev], *get_velocity_eb()[lev], 1);
         set_eb_density(lev, m_t_old[lev], *get_density_eb()[lev], 1);
         set_eb_tracer(lev, m_t_old[lev], *get_tracer_eb()[lev], 1);
       }
    }
#endif

    ApplyPredictor();

    if (m_advection_type == "MOL") {
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
        }
        amrex::Print() << "After Predictor - Befor Corrector" << "\n";
        ApplyCorrector();
    }
    // EY: Impose zero velocity for inside of pellet
    if (m_fluid_model == FluidModel::TwoMu)
    {
        // Update properties first 
        for(int lev = 0; lev <= finest_level; ++lev)
        {
           update_properties(lev, t_prop[lev]);
        }

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
                            {
                                // no internal flow
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
#ifdef INCFLO_USE_PARTICLES
    particleData.Redistribute();
#endif

#if 0
    // This sums over all levels
    if (m_test_tracer_conservation) {
        Real sum = volumeWeightedSum(get_tracer_new_const(),0,geom,ref_ratio);
        amrex::Print() << "Sum tracer volume wgt2 = " << m_cur_time+m_dt << " " <<
                           sum << std::endl;
    }
#endif

    // Stop timing current time step
    Real end_step = static_cast<Real>(ParallelDescriptor::second()) - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }
}

void incflo::update_properties(int lev, MultiFab& t_prop)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();

    int ncomp = t_prop.nComp();

    for (MFIter mfi(t_prop); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> t_prop_arr = t_prop.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            incflo::update_thermal_properties_and_phases(i, j, k, t_prop_arr);
        });
    }
}

