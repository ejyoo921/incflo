#ifndef _USERFUNCS_H_
#define _USERFUNCS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

using namespace amrex;
namespace sm_userfuncs
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void temperature_bc(int i, int j, int k,
                        int dir, int sgn,
                        Array4<Real> const& phi,
                        Array4<Real> const& bcarr,
                        Array4<Real> const& robin_a,
                        Array4<Real> const& robin_b,
                        Array4<Real> const& robin_f,
                        GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                        GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                        GpuArray<Real, AMREX_SPACEDIM> dx,
                        const Real time,
                        ProbParm const& prob_parm)
    {
        //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        const int im1 = (dir == 0) ? i-1 : i;
        const int jm1 = (dir == 1) ? j-1 : j;
        const int km1 = (dir == 2) ? k-1 : k;

        if(sgn == -1) 
        { // lo sides
            robin_a(im1,jm1,km1) = 1.0;
            robin_b(im1,jm1,km1) = 0.0;
            robin_f(im1,jm1,km1) = 0.0;
            bcarr(im1,jm1,km1) = 0.0;
        }
        else
        {
            robin_a(i,j,k) = 1.0;
            robin_b(i,j,k) = 0.0;
            robin_f(i,j,k) = 0.0;
            bcarr(i,j,k) = 0.0;
        } 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void update_thermal_properties_and_phases(int i, int j, int k,
                                   Array4<Real> const& phi,
                                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                   GpuArray<Real, AMREX_SPACEDIM> dx,
                                   const Real time,
                                   ProbParm const& prob_parm)
    {

        amrex::Real Temp=phi(i,j,k,TEMP_ID);
        amrex::Real cp_fe,cond_fe,dens_fe;
        amrex::Real cp_slg,cond_slg,dens_slg;
        amrex::Real sol_fe,mol_fe,sol_slg,mol_slg;

        //for iron
        if(Temp>prob_parm.Tmelt2_fe)
        {
            cp_fe     = prob_parm.cp_liq_fe;
            cond_fe   = prob_parm.k_liq_fe;
            dens_fe   = prob_parm.rho_liq_fe;
            sol_fe    = 0.0;
            mol_fe    = 1.0;
        }
        //solid
        else if(Temp<prob_parm.Tmelt1_fe)
        {
            cp_fe = prob_parm.cp_solid_fe;
            cond_fe = prob_parm.k_solid_fe;
            dens_fe = prob_parm.rho_solid_fe;
            sol_fe = 1.0;
            mol_fe = 0.0;
        }
        //phase change
        else
        {
            cp_fe = prob_parm.Cpmelt_fe;
            amrex::Real frac = (Temp-prob_parm.Tmelt1_fe)/(prob_parm.Tmelt2_fe-prob_parm.Tmelt1_fe);
            cond_fe = (1.0-frac)*prob_parm.k_solid_fe+frac*prob_parm.k_liq_fe;
            dens_fe = (1.0-frac)*prob_parm.rho_solid_fe+frac*prob_parm.rho_liq_fe;
            sol_fe = 1.0-frac;
            mol_fe = frac;
        }
        
        //for slag
        if(Temp>prob_parm.Tmelt2_slg)
        {
            cp_slg     = prob_parm.cp_liq_slg;
            cond_slg   = prob_parm.k_liq_slg;
            dens_slg   = prob_parm.rho_liq_slg;
            sol_slg    = 0.0;
            mol_slg    = 1.0;
        }
        //solid
        else if(Temp<prob_parm.Tmelt1_slg)
        {
            cp_slg   = prob_parm.cp_solid_slg;
            cond_slg = prob_parm.k_solid_slg;
            dens_slg = prob_parm.rho_solid_slg;
            sol_slg=1.0;
            mol_slg=0.0;
        }
        //phase change
        else
        {
            cp_slg = prob_parm.Cpmelt_slg;
            amrex::Real frac=(Temp-prob_parm.Tmelt1_slg)/(prob_parm.Tmelt2_slg-prob_parm.Tmelt1_slg);
            cond_slg  = (1.0-frac)*prob_parm.k_solid_slg + frac*prob_parm.k_liq_slg;
            dens_slg  = (1.0-frac)*prob_parm.rho_solid_slg + frac*prob_parm.rho_liq_slg;
            sol_slg = 1.0-frac;
            mol_slg = frac;
        }
        
        amrex::Real vfrac_fe=(phi(i,j,k,NTHERMVARS+MIXMASS_ID)-dens_slg)
        /(dens_fe-dens_slg);
        phi(i,j,k,DENS_ID)   = dens_slg*(1.0-vfrac_fe)+dens_fe*vfrac_fe;
        phi(i,j,k,COND_ID)   = cond_slg*(1.0-vfrac_fe)+cond_fe*vfrac_fe;
        phi(i,j,k,SPHEAT_ID)   = cp_slg*(1.0-vfrac_fe)+cp_fe*vfrac_fe;
        phi(i,j,k,NTHERMVARS+SOLFE_ID)   = sol_fe;             
        phi(i,j,k,NTHERMVARS+MOLFE_ID)   = mol_fe;             
        phi(i,j,k,NTHERMVARS+SOLSLG_ID)   = sol_slg;             
        phi(i,j,k,NTHERMVARS+MOLSLG_ID)   = mol_slg;             
        phi(i,j,k,NTHERMVARS+VFRAC_ID)  = vfrac_fe;

        if(phi(i,j,k,SPHEAT_ID)<0.0)
        {
            amrex::AllPrint()<<"sp heat less than 0:"<<vfrac_fe<<
            "\t"<<Temp<<"\t"<<phi(i,j,k,TEMP_ID)<<"\n";
            /*cp_slg<<"\t"<<cp_fe<<"\t"<<dens_fe<<"\t"<<dens_slg<<"\t"<<
            phi(i,j,k,MIXMASS_ID)<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\n";*/
            amrex::Abort();
        }
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void temperature_source(int i, int j, int k,
                            Array4<Real> const& phi,
                            Array4<Real> const& source,
                            GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                            GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                            GpuArray<Real, AMREX_SPACEDIM> dx,
                            const Real time,
                            ProbParm const& prob_parm)
    {
        source(i,j,k) += 0.0;
    }

}
#endif
