#include <incflo.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

using namespace amrex;

// tag cells for refinement
// overrides the pure virtual function in AmrCore
void incflo::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    BL_PROFILE("incflo::ErrorEst()");

    static bool first = true;
    static Vector<Real> rhoerr_v, gradrhoerr_v;

    // EY: tagging for Steelmake - jump viscosity
    static Vector<Real> etaerr_v, gradetaerr_v;
    // ----EY
    static bool tag_region;

    // only do this during the first call to ErrorEst
    if (first) {
        first = false;
        ParmParse pp("incflo");

        pp.queryarr("rhoerr", rhoerr_v);
        if (!rhoerr_v.empty()) {
            Real last = rhoerr_v.back();
            rhoerr_v.resize(max_level+1, last);
        }

        pp.queryarr("gradrhoerr", gradrhoerr_v);
        if (!gradrhoerr_v.empty()) {
            Real last = gradrhoerr_v.back();
            gradrhoerr_v.resize(max_level+1, last);
        }

        // EY: tagging for Steelmake - jump viscosity---
        pp.queryarr("etaerr", etaerr_v);
        if (!etaerr_v.empty()) {
            Real last = etaerr_v.back();
            etaerr_v.resize(max_level+1, last);
        }

        pp.queryarr("gradetaerr", gradetaerr_v);
        if (!gradetaerr_v.empty()) {
            Real last = gradetaerr_v.back();
            gradetaerr_v.resize(max_level+1, last);
        }

        // ----EY: More TODO
        tag_region_lo.resize(3);
        tag_region_hi.resize(3);

        tag_region = false;
        pp.query("tag_region", tag_region);

        pp.queryarr("tag_region_lo", tag_region_lo);
        pp.queryarr("tag_region_hi", tag_region_hi);
    }

    const auto   tagval = TagBox::SET;

    bool tag_rho = levc < rhoerr_v.size();
    bool tag_gradrho = levc < gradrhoerr_v.size();

    if (tag_gradrho) {
        fillpatch_density(levc, time, m_leveldata[levc]->density, 1);
    }

    // EY: tagging for Steelmake - jump viscosity
    bool tag_eta = levc < etaerr_v.size();
    bool tag_gradeta = levc < gradetaerr_v.size();
    // if (tag_gradeta) {
    //     fillpatch_viscosity(lev, time, m_leveldata[lev]->viscosity, 1);
    // }

    AMREX_D_TERM(const Real l_dx = geom[levc].CellSize(0);,
                 const Real l_dy = geom[levc].CellSize(1);,
                 const Real l_dz = geom[levc].CellSize(2););

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_leveldata[levc]->density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);
        
        if (tag_rho || tag_gradrho) 
        {
            Array4<Real const> const& rho = m_leveldata[levc]->density.const_array(mfi);
            Real rhoerr = tag_rho ? rhoerr_v[levc]: std::numeric_limits<Real>::max();
            Real gradrhoerr = tag_gradrho ? gradrhoerr_v[levc] : std::numeric_limits<Real>::max();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (tag_rho && rho(i,j,k) > rhoerr) {
                    tag(i,j,k) = tagval;
                }
                if (tag_gradrho) {
                    Real ax = amrex::Math::abs(rho(i+1,j,k) - rho(i,j,k));
                    Real ay = amrex::Math::abs(rho(i,j+1,k) - rho(i,j,k));
                    ax = amrex::max(ax,amrex::Math::abs(rho(i,j,k) - rho(i-1,j,k)));
                    ay = amrex::max(ay,amrex::Math::abs(rho(i,j,k) - rho(i,j-1,k)));
#if (AMREX_SPACEDIM == 2)
                    if (amrex::max(ax,ay) >= gradrhoerr) {
                        tag(i,j,k) = tagval;
                    }
#elif (AMREX_SPACEDIM == 3)
                    Real az = amrex::Math::abs(rho(i,j,k+1) - rho(i,j,k));
                    az = amrex::max(az,amrex::Math::abs(rho(i,j,k) - rho(i,j,k-1)));
                    if (amrex::max(ax,ay,az) >= gradrhoerr) {
                        tag(i,j,k) = tagval;
                    }
#endif
                }
            });
        }

        if (tag_region) {

            Real xlo = tag_region_lo[0];
            Real ylo = tag_region_lo[1];
            Real xhi = tag_region_hi[0];
            Real yhi = tag_region_hi[1];
            auto const& problo = geom[levc].ProbLoArray();

#if (AMREX_SPACEDIM == 2)

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 Real x = problo[0] + (i+0.5)*l_dx;
                 Real y = problo[1] + (j+0.5)*l_dy;

                 // Tag if we are inside the specified box
                 if (x >= xlo && x <= xhi && y >= ylo && y <= yhi)
                 {
                    tag(i,j,k) = tagval;
                 }
            });

#else
            Real zlo = tag_region_lo[2];
            Real zhi = tag_region_hi[2];

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 Real x = problo[0] + Real(i+0.5)*l_dx;
                 Real y = problo[1] + Real(j+0.5)*l_dy;
                 Real z = problo[2] + Real(k+0.5)*l_dz;

                 // Tag if we are inside the specified box
                 if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
                 {
                    tag(i,j,k) = tagval;
                 }
            });
#endif
        }
    } // mfi

    // // EY: Tagging Eta (viscosity)
    for (MFIter mfi(m_leveldata[lev]->viscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);

        // EY: we may need one more tag here for Steelmaking (|| tag_)
        if (tag_eta || tag_gradeta) 
        {
            // //EY: 
            Array4<Real const> eta = m_leveldata[lev]->viscosity.const_array(mfi);
            Real etaerr = tag_eta ? etaerr_v[lev]: std::numeric_limits<Real>::max();
            Real gradetaerr = tag_gradeta ? gradetaerr_v[lev] : std::numeric_limits<Real>::max();
            // EY: Need to implement leveldata for Eta (viscosity) -> done
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (tag_eta && eta(i,j,k) > etaerr) {
                    tag(i,j,k) = tagval;
                }
                if (tag_gradeta) {
                    Real ax = amrex::Math::abs(eta(i+1,j,k) - eta(i,j,k));
                    Real ay = amrex::Math::abs(eta(i,j+1,k) - eta(i,j,k));
                    ax = amrex::max(ax,amrex::Math::abs(eta(i,j,k) - eta(i-1,j,k)));
                    ay = amrex::max(ay,amrex::Math::abs(eta(i,j,k) - eta(i,j-1,k)));
#if (AMREX_SPACEDIM == 2)
                    if (amrex::max(ax,ay) >= gradetaerr) {
                        tag(i,j,k) = tagval;
                    }
#elif (AMREX_SPACEDIM == 3)
                    Real az = amrex::Math::abs(eta(i,j,k+1) - eta(i,j,k));
                    az = amrex::max(az,amrex::Math::abs(eta(i,j,k) - eta(i,j,k-1)));
                    
                    if (amrex::max(ax,ay,az) >= gradetaerr)
                    {
                        tag(i,j,k) = tagval;
                    }
#endif
                }
            });
        }

        if (tag_region) {

            Real xlo = tag_region_lo[0];
            Real ylo = tag_region_lo[1];
            Real xhi = tag_region_hi[0];
            Real yhi = tag_region_hi[1];
            auto const& problo = geom[lev].ProbLoArray();

#if (AMREX_SPACEDIM == 2)

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 Real x = problo[0] + (i+0.5)*l_dx;
                 Real y = problo[1] + (j+0.5)*l_dy;

                 // Tag if we are inside the specified box
                 if (x >= xlo && x <= xhi && y >= ylo && y <= yhi)
                 {
                    tag(i,j,k) = tagval;
                 }
            });

#else
            Real zlo = tag_region_lo[2];
            Real zhi = tag_region_hi[2];

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 Real x = problo[0] + Real(i+0.5)*l_dx;
                 Real y = problo[1] + Real(j+0.5)*l_dy;
                 Real z = problo[2] + Real(k+0.5)*l_dz;

                 // Tag if we are inside the specified box
                 if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
                 {
                    tag(i,j,k) = tagval;
                 }
            });
#endif
        }
    }

#ifdef AMREX_USE_EB
    m_refine_cutcells = true;
    // Refine on cut cells
    if (m_refine_cutcells)
    {
        amrex::TagCutCells(tags, m_leveldata[levc]->velocity);
    }
#endif

#ifdef INCFLO_USE_PARTICLES
    if (m_refine_particles)
    {
        //
        // This allows dynamic refinement based on the number of particles per cell
        //
        // Note that we must count all the particles in levels both at and above the current,
        //      since otherwise, e.g., if the particles are all at level 1, counting particles at
        //      level 0 will not trigger refinement when regridding so level 1 will disappear,
        //      then come back at the next regridding
        //
        const auto& particles_namelist( particleData.getNames() );
        std::unique_ptr<MultiFab> mf = std::make_unique<MultiFab>(grids[levc], dmap[levc], 1, 0);
        mf->setVal(0.0);
        IntVect rr = IntVect::TheUnitVector();
        for (int lev = levc; lev <= finest_level; lev++)
        {
            MultiFab temp_dat(grids[lev], dmap[lev], 1, 0); temp_dat.setVal(0);
            particleData[particles_namelist[0]]->IncrementWithTotal(temp_dat, lev);

            MultiFab temp_dat_crse(grids[levc], dmap[levc], 1, 0); temp_dat_crse.setVal(0);

            if (lev == levc) {
                MultiFab::Copy(*mf, temp_dat, 0, 0, 1, 0);
            } else {
                for (int d = 0; d < AMREX_SPACEDIM; d++) {
                    rr[d] *= ref_ratio[levc][d];
                }
                average_down(temp_dat, temp_dat_crse, 0, 1, rr);
                MultiFab::Add(*mf, temp_dat_crse, 0, 0, 1, 0);
            }
        } // lev
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_leveldata[levc]->density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            auto const&  mf_arr = mf->const_array(mfi);
            auto const& tag_arr = tags.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mf_arr(i,j,k) > 0) {
                    tag_arr(i,j,k) = tagval;
                }
            });
        } // mfi
    } // if m_refine_particles
#endif
}
