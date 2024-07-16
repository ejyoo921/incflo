#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <ProbParm.H>
#include <Constants.H>

using namespace amrex;

AMREX_INLINE
void amrex_probinit (ProbParm *h_prob_parm, ProbParm *d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");

    // Extract position and velocities
    amrex::Vector<amrex::Real> rads;
    amrex::Vector<amrex::Real> centx;
    amrex::Vector<amrex::Real> centy;
    amrex::Vector<amrex::Real> centz;
    
    pp.get("npellets", h_prob_parm->npellets );
    
    rads.resize(h_prob_parm->npellets);
    centx.resize(h_prob_parm->npellets);
    centy.resize(h_prob_parm->npellets);
    centz.resize(h_prob_parm->npellets);

    pp.getarr("pellet_rads",  rads);    
    pp.getarr("pellet_centx", centx);
    pp.getarr("pellet_centy", centy);
    pp.getarr("pellet_centz", centz);

    pp.query("Tinit_solid",h_prob_parm->Tinit_solid_fe);
    pp.query("Tinit_liq",h_prob_parm->Tinit_liq_slg);

    h_prob_parm->pellet_rads  = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm->npellets*sizeof(amrex::Real));
    h_prob_parm->pellet_centx = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm->npellets*sizeof(amrex::Real));
    h_prob_parm->pellet_centy = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm->npellets*sizeof(amrex::Real));
    h_prob_parm->pellet_centz = (amrex::Real*) amrex::The_Arena()->alloc(h_prob_parm->npellets*sizeof(amrex::Real));

    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                rads.begin(),
                                rads.end(),
                       h_prob_parm->pellet_rads);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                centx.begin(),
                                centx.end(),
                       h_prob_parm->pellet_centx);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                centy.begin(),
                                centy.end(),
                       h_prob_parm->pellet_centy);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                centz.begin(),
                                centz.end(),
                       h_prob_parm->pellet_centz);

    amrex::Gpu::copy(Gpu::hostToDevice, h_prob_parm, h_prob_parm+1,d_prob_parm);

}

AMREX_GPU_DEVICE
AMREX_INLINE
void initData(Box const& bx, Array4<Real> const& phi, 
              GeometryData const& geomdata,
              ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            for (int i = lo.x; i <= hi.x; ++i) 
            {

                Real vfrac_fe=0.0;

                for(int kk=0;kk<2;kk++)
                {
                    for(int jj=0;jj<2;jj++)
                    {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real x = prob_lo[0] + (i+ii) * dx[0];
                            Real y = prob_lo[1] + (j+jj) * dx[1];
                            Real z = prob_lo[2] + (k+kk) * dx[2];

                            for(int np=0;np<prob_parm->npellets;np++)
                            {

                                Real dist2 = std::pow(x-prob_parm->pellet_centx[np],2.0)+                
                                std::pow(y-prob_parm->pellet_centy[np],2.0)+                
                                std::pow(z-prob_parm->pellet_centz[np],2.0); 

                                if(dist2 < std::pow(prob_parm->pellet_rads[np],2.0))
                                {              
                                    vfrac_fe+=1.0;
                                    break;
                                }
                            }
                        }
                    }
                }

                vfrac_fe=vfrac_fe/8.0;

                Real sol_fe,mol_fe,sol_slg,mol_slg;

                sol_fe=vfrac_fe;
                mol_fe=0.0;
                sol_slg=0.0; 
                mol_slg=1.0-vfrac_fe;

                phi(i,j,k,TEMP_ID)  =  prob_parm->Tinit_solid_fe*sol_fe
                + prob_parm->Tinit_liq_slg*mol_slg;

                phi(i,j,k,DENS_ID)  = prob_parm->rho_solid_fe*sol_fe
                + prob_parm->rho_liq_slg*mol_slg;

                phi(i,j,k,SPHEAT_ID)  =   prob_parm->cp_solid_fe*sol_fe
                + prob_parm->cp_liq_slg*mol_slg;

                phi(i,j,k,COND_ID)  = prob_parm->k_solid_fe*sol_fe +
                prob_parm->k_liq_slg*mol_slg;

                phi(i,j,k,NTHERMVARS+SOLFE_ID)    = sol_fe;             
                phi(i,j,k,NTHERMVARS+MOLFE_ID)    = mol_fe;             
                phi(i,j,k,NTHERMVARS+SOLSLG_ID)   = sol_slg;             
                phi(i,j,k,NTHERMVARS+MOLSLG_ID)   = mol_slg;             
                phi(i,j,k,NTHERMVARS+MIXMASS_ID)   = phi(i,j,k,DENS_ID);     
                if(phi(i,j,k,NTHERMVARS+MIXMASS_ID)<1.0)
                {
                    amrex::Print()<<"mixmass:"<<phi(i,j,k,NTHERMVARS+MIXMASS_ID)<<"\n";
                }        
                phi(i,j,k,NTHERMVARS+VFRAC_ID)   = vfrac_fe;             
            }
        }
    }
}
#endif
