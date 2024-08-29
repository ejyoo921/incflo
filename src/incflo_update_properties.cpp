#include <incflo.H>

using namespace amrex;

namespace 
{
    // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    // amrex::Real compute_vfrac (int lev, int i, int j, int k, 
    //                             amrex::GpuArray<double, 3> prob_lo, 
    //                             amrex::GpuArray<double, 3> prob_hi, 
    //                             amrex::GpuArray<double, 3> dx)
    // {
    //     amrex::ParmParse pp("prob");
    //     // Extract position and velocities
    //     int npellets = 0;
    //     amrex::Vector<amrex::Real> rads;
    //     amrex::Vector<amrex::Real> centx;
    //     amrex::Vector<amrex::Real> centy;
    //     amrex::Vector<amrex::Real> centz;

    //     pp.get("npellets", npellets);
    //     pp.getarr("pellet_rads",  rads);    
    //     pp.getarr("pellet_centx", centx);
    //     pp.getarr("pellet_centy", centy);
    //     pp.getarr("pellet_centz", centz);

    // #ifdef _OPENMP
    // #pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
    // #endif

    //     Real vfrac_fe = 0.0;

    //     for(int kk=0;kk<2;kk++)
    //     {
    //         for(int jj=0;jj<2;jj++)
    //         {
    //             for(int ii=0;ii<2;ii++)
    //             {
    //                 Real x = prob_lo[0] + (i+ii) * dx[0];
    //                 Real y = prob_lo[1] + (j+jj) * dx[1];
    //                 Real z = prob_lo[2] + (k+kk) * dx[2];

    //                 for(int np = 0; np < npellets; np++)
    //                 {
    //                     Real dist2 = std::pow(x - centx[np], 2.0)+                
    //                                 std::pow(y - centy[np], 2.0)+                
    //                                 std::pow(z - centz[np], 2.0); 

    //                     if(dist2 < std::pow(rads[np], 2.0))
    //                     {
    //                         amrex::Print() << "vfrac add up?" << "\n";          
    //                         vfrac_fe+=1.0;
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     vfrac_fe=vfrac_fe/8.0;

    //     return vfrac_fe;        
    // }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real bound01(const amrex::Real q0)
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
    amrex::Real visc_0(const amrex::Real Temp)
    { // artificial viscosity for simultation
        amrex::Real visc_1a = 0.01;
        amrex::Real MeltTemp = 800.0;
        if (Temp >= MeltTemp)
        { 
            visc_1a = 0.01;
        }
        else
        { 
            visc_1a = 100.0;
        }
        return visc_1a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real cp_0(const amrex::Real Temp)
    { //cp calc for solid pellet phase 0, a test material used to exercise the code  
     
        amrex::Real cp_0a = 0.400;
        amrex::Real cp_1  = 0.400;
        amrex::Real cp_2  = 0.600;
        amrex::Real T1    = 800.0;
        amrex::Real TL    = 790.0;
        amrex::Real TH    = 810.0;
        amrex::Real S     = (cp_2-cp_1)/(TH-TL);
        amrex::Real omega = 10.0;
        amrex::Real c     = 10.0;
      
     
        if (Temp <= TL)
        { 
            cp_0a = cp_1+ c * exp(-0.5*pow(((Temp-T1)/omega),2));
        }
        else
        { if (Temp <= TH)
            { 
                cp_0a = cp_1 + S*(Temp-TL) + c * exp(-0.5*pow(((Temp-T1)/omega),2));
            }
            else
            { 
                cp_0a = cp_2 +c * exp(-0.5*pow(((Temp-T1)/omega),2));
            }
        }

        return 1000.0*cp_0a;
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real k_3(const amrex::Real Temp)
    {   //k calc for 1 % C bath
        amrex::Real k3a = 30.0;
        amrex::Real MeltTemp = 1669.0;
        if (Temp >= MeltTemp)
        { k3a = 30.0;
        }
        else
        { k3a =39.0;
        }
        return k3a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real k_4(const amrex::Real Temp) 
    {   //k calc for 3% carbon bath 
        amrex::Real k_4a = 39.0;
        amrex::Real MeltTemp = 1477.0;
        if (Temp >= MeltTemp)
        { k_4a =30.0;
        }
        else
        { k_4a = 39.0;
        }
        return k_4a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real k_5(const amrex::Real Temp)
    {   //k calc for slag  
        amrex::Real k_5a = 1.0;
        amrex::Real MeltTemp = 1574.0;
        if (Temp >= MeltTemp)
        { k_5a = 0.1;
        }
        else
        { k_5a = 1.0;
        }
        return k_5a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rho_3(const amrex::Real Temp)
    {   //rho calc for 1 % C bath
        amrex::Real rho3a = 7800;
        amrex::Real MeltTemp = 1669.0;
        if (Temp >= MeltTemp)
        { rho3a = 6900.0;
        }
        else
        { rho3a =7800.0;
        }
            return rho3a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rho_4(const amrex::Real Temp)
    {   //rho calc for 3% carbon bath  
        amrex::Real k_4a = 39.0;
        amrex::Real MeltTemp = 1477.0;
        if (Temp >= MeltTemp)
        {   k_4a =30.0;
        }
        else
        {   k_4a = 39.0;
        }
        return k_4a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rho_5(const Real Temp)
    {   //rho calc for slag  
        amrex::Real rho_5a = 2700.0;
        amrex::Real MeltTemp = 1574.0;
        if (Temp >= MeltTemp)
        {   rho_5a = 2700.0;
        }
        else
        {   rho_5a = 3000.0;
        }
        return rho_5a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_0(const amrex::Real Temp)
    { // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_0a = 0.0;
        amrex::Real MeltTemp = 800.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_0a = 1.0;
        }
        else
        { 
            liqfrac_0a = 0.0;
        }
        return liqfrac_0a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_1(const amrex::Real Temp)
    { // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_1a = 0.0;
        amrex::Real MeltTemp = 1806.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_1a = 1.0;
        }
        else
        { 
            liqfrac_1a = 0.0;
        }
        return liqfrac_1a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_2(const amrex::Real Temp)
    {  // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_2a = 0.0;
        amrex::Real MeltTemp = 1427.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_2a = 1.0;
        }
        else
        { 
            liqfrac_2a = 0.0;
        }
        return liqfrac_2a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_3(const amrex::Real Temp)
    { // calculation for 1% carbon bath
      // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_3a = 0.0;
        amrex::Real MeltTemp = 1669.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_3a = 1.0;
        }
        else
        { 
            liqfrac_3a = 0.0;
        }
        return liqfrac_3a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_4(const amrex::Real Temp)
    { // calculation for 3% carbon bath
      // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_4a = 0.0;
        amrex::Real MeltTemp = 1477.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_4a = 1.0;
        }
        else
        { 
            liqfrac_4a = 0.0;
        }
        return liqfrac_4a;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real liqfrac_5(const amrex::Real Temp)
    { // phase 5 values are for slag  
      // a value of 1 indicates liquid, a value of 0 is solid.
        amrex::Real liqfrac_5a = 0.0;
        amrex::Real MeltTemp = 1574.0;
        if (Temp >= MeltTemp)
        { 
            liqfrac_5a = 1.0;
        }
        else
        { 
            liqfrac_5a = 0.0;
        }
        return liqfrac_5a;
    }

    //******************* Get properties here ******************
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real compute_cp(const amrex::Real Temp, int phase)
    {   // calculate cp based on phase and temp
        // unit: J/K/kg
        Real cpcalc = 1000.0;
        switch (phase) {
            case 0:
            cpcalc = cp_0(Temp);
            break;

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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real compute_k(const amrex::Real Temp,int phase)
    { // calculate rho based on phase and temp
        Real kcalc = 1000.0; 
        switch (phase) {
            case 0:
            kcalc = 36.5;
	        break;

            case 1:
            kcalc = 36.5;
            break;

            case 2:
            kcalc = 2.90;
            break;

            case 3:
            kcalc = k_3(Temp);
            break;

            case 4:
            kcalc = k_4(Temp);
            break;

            case 5:
            kcalc = k_5(Temp);
            break;
        }
        return kcalc;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real compute_rho(const amrex::Real Temp,int phase)
    {   // calculate rho based on phase and temp
        Real rhocalc=1000.0;
        switch (phase) {
            case 0:
                rhocalc = 7900.0;
            break;

            case 1:
                rhocalc = 7900.0;
            break;

            case 2:
                rhocalc = 7700.0;
            break;

            case 3:
                rhocalc = rho_3(Temp);
            break;

            case 4:
                rhocalc = rho_4(Temp);
            break;

            case 5:
                rhocalc = rho_5(Temp);
            break;
        }
        return rhocalc;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real compute_liqfrac(const amrex::Real Temp,int phase)
    { // calculate liqfrac based on phase and temp
        Real liqfraccalc=0.0;
        int material1 = phase;
        switch (material1) {
	    case 0:
           liqfraccalc = liqfrac_0(Temp);
        break;

        case 1:
           liqfraccalc = liqfrac_1(Temp);
        break;

        case 2:
           liqfraccalc = liqfrac_2(Temp);
        break;

        case 3:
           liqfraccalc = liqfrac_3(Temp);
        break;

        case 4:
           liqfraccalc = liqfrac_4(Temp);
        break;

        case 5:
           liqfraccalc = liqfrac_5(Temp);
        break;
      }
      return liqfraccalc;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_visc(const amrex::Real Temp,int phase)
    {  // calculate viscosity based on phase and temp
        amrex::Real visccalc=1000.0;
        switch (phase) {
        case 0:
            visccalc = visc_0(Temp);
        break;
        case 1:
            visccalc = 0.02;
        break;
        case 2:
            visccalc = 0.02;
        break;
        case 3:
            visccalc = 0.02;
        break;
        case 4:
            visccalc = 0.02;
        break;
        case 5:
            visccalc = 0.02;
        break;
      }
      return visccalc;
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
        
        auto& ld = *m_leveldata[lev];
        for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& gbx = mfi.growntilebox(); //bigger box (with ghosts)

            Array4<Real> vfrac_mix_arr = ld.vfrac_mix.array(mfi); 
            // Array4<Real> temp_arr = ld.tracer.array(mfi);
            Array4<Real> cp_arr   = ld.cp_steel.array(mfi); 
            Array4<Real> dens_arr = ld.rho_steel.array(mfi); 
            Array4<Real> cond_arr = ld.k_steel.array(mfi); 
            Array4<Real> eta_arr = ld.viscosity.array(mfi); 
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& temp_arr   = ld.tracer_o.const_array(mfi);

            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // incflo::update_thermal_properties_and_phases(i, j, k, t_prop_arr, prob_lo, prob_hi, dx);
                // For SteelSlag case 
                amrex::Real Temp = temp_arr(i,j,k); // Need to bring actual Temp vals
                amrex::Real vfrac_fe_mix = vfrac_mix_arr(i,j,k);
                amrex::Real vfrac_fe;
                amrex::Real cp_fe,  cond_fe,  dens_fe; // Fe
                amrex::Real cp_slg, cond_slg, dens_slg; // Slag
                amrex::Real sol_fe, mol_fe, sol_slg, mol_slg; // Phases

                const auto lo = lbound(gbx);
                const auto hi = ubound(gbx);    

                auto const& prob_lo = geom[lev].ProbLoArray();
                auto const& prob_hi = geom[lev].ProbLoArray();
                auto const& dx = geom[lev].CellSizeArray();  

                // TODO: This should not be computed for every time step.
                // Just once when initializing and keep it somewhere.
                // vfrac_fe_mix = compute_vfrac (lev, i, j, k, prob_lo, prob_hi, dx);

                // update rho ----------------------------------------------------
                dens_fe         = compute_rho(Temp,0); //zero is temporary
                dens_slg        = compute_rho(Temp,5);


                vfrac_fe        = bound01((vfrac_fe_mix-dens_slg)/(dens_fe-dens_slg));
                dens_arr(i,j,k)  = dens_slg*(1.0-vfrac_fe) + dens_fe*vfrac_fe;
                // dens_arr(i,j,k)  = 1.0;

                // update cp -----------------------------------------------------
                cp_fe           = compute_cp(Temp, 0); //zero is temporary
                cp_slg          = compute_cp(Temp, 5);
                cp_arr(i,j,k)   = cp_slg*(1.0-vfrac_fe) + cp_fe*vfrac_fe;
                cp_arr(i,j,k)  = 1.0;

                // update conductivity -------------------------------------------
                cond_fe         = compute_k(Temp,0); //zero is temporary
                cond_slg        = compute_k(Temp,5); //k is conductivity
                cond_arr(i,j,k) = cond_slg*(1.0-vfrac_fe) + cond_fe*vfrac_fe;
                cond_arr(i,j,k)  = 0.0;

                // get iron properties
	            mol_fe = bound01(compute_liqfrac(Temp,0));
                sol_fe = bound01((1.0 - mol_fe));

                if (i == 7 & j == 7 & k == 7)
                {
                    amrex::Print() << "Temperature = " << Temp << "\n";
                    amrex::Print() << "VFRAC = " << vfrac_fe << "\n";
                    amrex::Print() << "Cp = " << cp_arr(i,j,k) << "\n";
                    amrex::Print() << "Conductivity = " << cond_arr(i,j,k) << "\n";
                    amrex::Print() << "Density = " << dens_arr(i,j,k) << "\n";

                }

                // update phases -------------------------------------------------
                // phi(i,j,k,NTHERMVARS+SOLFE_ID)   = vfrac_fe*sol_fe;             
                // phi(i,j,k,NTHERMVARS+MOLFE_ID)   = vfrac_fe*mol_fe;             
                // phi(i,j,k,NTHERMVARS+SOLSLG_ID)  = bound01((1.0-vfrac_fe))*sol_slg;             
                // phi(i,j,k,NTHERMVARS+MOLSLG_ID)  = bound01((1.0-vfrac_fe))*mol_slg;  

                // set zero velocity - inside pellet
                ParmParse pp("incflo");
                bool m_zero_vel = false;
                pp.queryAdd("zero_vel", m_zero_vel);

                if (m_zero_vel) // Make zero velocity 
                {
                    int npellets = 0;
                    amrex::ParmParse pp("prob");
                    pp.get("npellets", npellets);

                    int inside_pellet = 0;
                    for(int np = 0; np < npellets; np++)
                    {
                        // Update the inside location with Temperature
                        Real Temp_melt = 1500.0;
                        if (Temp < Temp_melt)
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
                        eta_arr(i,j,k) = m_mu*pow(10, m_n_0);
                    } // inside pellet
                } // if-zero-vel
            }); // i,j,k
        } // mfi
    } // lev

}
