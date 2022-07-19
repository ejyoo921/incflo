#include <incflo.H>
#include <incflo_derive_K.H>

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expterm (amrex::Real nu) noexcept
{
    return (nu < 1.e-9) ? (1.0-0.5*nu+nu*nu*(1.0/6.0)-(nu*nu*nu)*(1./24.))
                        : -std::expm1(-nu)/nu; 
}

// Compute the I term, where I is the inertial number, in mu(I) relation
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertialNum (amrex::Real sr, amrex::Real p_ext, amrex::Real ro_0, amrex::Real diam, amrex::Real mu_1, amrex::Real A_1, amrex::Real alpha_1) noexcept
{
    return mu_1 + A_1*std::pow((sr/2)*diam/std::pow(p_ext/ro_0, 0.5), alpha_1);
}

struct NonNewtonianViscosity //Apparent viscosity
{
    incflo::FluidModel fluid_model;
    amrex::Real mu, n_flow, tau_0, eta_0, papa_reg, ro_0, p_bg, diam, mu_1, A_1, alpha_1, mu_2, A_2, alpha_2;
    

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real sr, amrex::Real p_ext)const noexcept {
        switch (fluid_model)
        {
        case incflo::FluidModel::powerlaw:
        {
            return mu * std::pow(sr,n_flow-1.0);
        }
        case incflo::FluidModel::Bingham:
        {   
            amrex::Print() << "are we here"  << "\n";
            return mu + tau_0 * expterm(sr/papa_reg) / papa_reg;
            // expterm(sr/papa_reg)/papa_reg = (1-exp(-sr/papa_reg))/sr --> 1/2dot(gamma)'s papa ver.
        }
        case incflo::FluidModel::HerschelBulkley:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr/papa_reg)/papa_reg;
        }
        case incflo::FluidModel::deSouzaMendesDutra:
        {
            return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr*(eta_0/tau_0))*(eta_0/tau_0);
        }
        case incflo::FluidModel::Granular:
        {
            // amrex::Print() << "mu_1 = " << mu_1 << "\n";
            // If you want pressure gradient add p_ext to p_bg
            return 2*(expterm(sr/papa_reg) / papa_reg)*(p_bg+p_ext)*inertialNum(sr, p_bg+p_ext, ro_0, diam, mu_1, A_1, alpha_1);
            //return std::pow(2*(expterm(sr/papa_reg) / papa_reg),2)*(p_bg+p_ext)*inertialNum(sr, p_bg+p_ext, ro_0, diam, mu_2, A_2, 2*alpha_2);
        }
        case incflo::FluidModel::Granular2:
        {
            // amrex::Print() << "mu_2 = "<< mu_2 <<"\n";
            return std::pow(2*(expterm(sr/papa_reg) / papa_reg),2)*(p_bg+p_ext)*inertialNum(sr, p_bg+p_ext, ro_0, diam, mu_2, A_2, 2*alpha_2);
        }
        default:
        {
            return mu;
        }
        };
    }
};

}

void incflo::compute_viscosity (Vector<MultiFab*> const& vel_eta,
                                Vector<MultiFab*> const& rho,
                                Vector<MultiFab*> const& vel,
                                Real time, int nghost)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        compute_viscosity_at_level(lev, vel_eta[lev], rho[lev], vel[lev], geom[lev], time, nghost);
    }
}

void incflo::compute_viscosity_at_level (int lev,
                                         MultiFab* vel_eta,
                                         MultiFab* /*rho*/,
                                         MultiFab* vel,
                                         Geometry& lev_geom,
                                         Real /*time*/, int nghost)
{
    if (m_fluid_model == FluidModel::Newtonian)
    {
        vel_eta->setVal(m_mu, 0, 1, nghost);
    }
    else if (m_fluid_model == FluidModel::Granular || m_fluid_model == FluidModel::Granular2)
    {
        /* code */
        NonNewtonianViscosity non_newtonian_viscosity;
        non_newtonian_viscosity.fluid_model = m_fluid_model;
        non_newtonian_viscosity.mu = m_mu;
        non_newtonian_viscosity.n_flow = m_n_0;
        non_newtonian_viscosity.tau_0 = m_tau_0;
        non_newtonian_viscosity.eta_0 = m_eta_0;
        non_newtonian_viscosity.papa_reg = m_papa_reg;

        //EY: for granular rheology
        non_newtonian_viscosity.ro_0 = m_ro_0;
        non_newtonian_viscosity.diam = m_diam;

        non_newtonian_viscosity.mu_1 = m_mu_1;
        non_newtonian_viscosity.A_1 = m_A_1;
        non_newtonian_viscosity.alpha_1 =m_alpha_1;

        non_newtonian_viscosity.mu_2 = m_mu_2;
        non_newtonian_viscosity.A_2 = m_A_2;
        non_newtonian_viscosity.alpha_2 =m_alpha_2;

        non_newtonian_viscosity.p_bg =m_p_bg;


        Real idx = 1.0 / lev_geom.CellSize(0);
        Real idy = 1.0 / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
        Real idz = 1.0 / lev_geom.CellSize(2);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta->array(mfi);
                Array4<Real const> const& vel_arr = vel->const_array(mfi);

                {
                    //EY: 1st attempt
                    amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);

                        Real nn = (m_vert_hi - m_vert_lo)/m_vert_n;
                        Real p_ext = m_gp0[1]*(j*nn);

                        eta_arr(i,j,k) = non_newtonian_viscosity(sr, p_ext);
                    });
        
                }
        }
    }
    else
    {
        NonNewtonianViscosity non_newtonian_viscosity;
        non_newtonian_viscosity.fluid_model = m_fluid_model;
        non_newtonian_viscosity.mu = m_mu;
        non_newtonian_viscosity.n_flow = m_n_0;
        non_newtonian_viscosity.tau_0 = m_tau_0;
        non_newtonian_viscosity.eta_0 = m_eta_0;
        non_newtonian_viscosity.papa_reg = m_papa_reg;

#ifdef AMREX_USE_EB
        auto const& fact = EBFactory(lev);
        auto const& flags = fact.getMultiEBCellFlagFab();
#endif

        Real idx = 1.0 / lev_geom.CellSize(0);
        Real idy = 1.0 / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
        Real idz = 1.0 / lev_geom.CellSize(2);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta->array(mfi);
                Array4<Real const> const& vel_arr = vel->const_array(mfi);
#ifdef AMREX_USE_EB
                auto const& flag_fab = flags[mfi];
                auto typ = flag_fab.getType(bx);
                if (typ == FabType::covered)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        eta_arr(i,j,k) = 0.0;
                    });
                }
                else if (typ == FabType::singlevalued)
                {
                    auto const& flag_arr = flag_fab.const_array();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate_eb(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr,flag_arr(i,j,k));

                        Real p_ext = 0.;
                        eta_arr(i,j,k) = non_newtonian_viscosity(sr, p_ext);
                    });
                }
                else
#endif
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        Real sr = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);
                        Real p_ext = 0.;

                        eta_arr(i,j,k) = non_newtonian_viscosity(sr, p_ext);
                    });
                }
        }
    }
}

void incflo::compute_tracer_diff_coeff (Vector<MultiFab*> const& tra_eta, int nghost)
{
    for (auto mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);
        }
    }
}
