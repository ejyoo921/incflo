#ifndef _USERFUNCS_H_
#define _USERFUNCS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

using namespace amrex;
namespace sm_userfuncs
{


  //  Start functions to calculate properties
    
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real bound01(const Real q0)
    //bounds a quantity between 0 and 1
    {
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


    //  Calculation of properties for phase 0, a test phase

      AMREX_GPU_DEVICE AMREX_INLINE

    Real cp_0(const Real Temp1)
    //cp calc for solid pellet phase 0, a test material used to exercise the code  
    {
     amrex::Real Temp = Temp1 + 273.15;
      amrex::Real cp_1a = 0.400;
      amrex::Real cp_1  = 0.400;
      amrex::Real cp_2  = 0.600;
      amrex::Real T1    = 800.0;
      amrex::Real TL    = 790.0;
      amrex::Real TH    = 810.0;
      amrex::Real S     = (cp_2-cp_1)/(TH-TL);
      amrex::Real omega = 10.0;
      amrex::Real c     = 10.0;
      amrex::Real Temp = Temp1 + 273.15;
      
     
      if (Temp <= TL)
      	{ cp_1a = cp_1+ c * exp(-0.5*pow(((Temp-T1)/omega),2));
      	}
      else
        { if (Temp <= TH)
            { cp_1a = cp_1 + S*(Temp-TL) + c * exp(-0.5*pow(((Temp-T1)/omega),2));
            }
          else
      	    { cp_1a = cp_2 +c * exp(-0.5*pow(((Temp-T1)/omega),2));
	    }}
       
      return 1000.0*cp_1a;
    }
  
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real k_0(const Real Temp1)
    //k calc for test phase 0  
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 36.5;
    }
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_0(const Real Temp1)
    //rho calc for solid pellet phase 1  
    {
     amrex::Real Temp = Temp1 + 273.15;
      return 7900.0;
    }

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_0(const Real Temp1)
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_1a = 0.0;
      amrex::Real MeltTemp = 800.0;
      if (Temp >= MeltTemp)
	{ liqfrac_1a = 1.0;
	}
      else
	{ liqfrac_1a = 0.0;
	}
      return liqfrac_1a;
    }


      AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_0(const Real Temp1)
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real visc_1a = 0.01;
      amrex::Real MeltTemp = 800.0;
      if (Temp >= MeltTemp)
	{ visc_1a = 0.01;
	}
      else
	{ visc_1a = 100.0;
	}
      return visc_1a;
    }


    //  Calculation of properties for phase 1
  
    AMREX_GPU_DEVICE AMREX_INLINE

     Real cp_1(const Real Temp1)
    //cp calc for solid pellet phase 1  
    {
      amrex::Real Temp = Temp1 + 273.15;
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
    
    Real k_1(const Real Temp1)
    //k calc for solid pellet phase 1  
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 36.5;
    }
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_1(const Real Temp1)
    //rho calc for solid pellet phase 1  
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 7900.0;
    }

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_1(const Real Temp1)
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_1a = 0.0;
      amrex::Real MeltTemp = 1806.0;
      if (Temp >= MeltTemp)
	{ liqfrac_1a = 1.0;
	}
      else
	{ liqfrac_1a = 0.0;
	}
      return liqfrac_1a;
    }

      AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_1(const Real Temp1)
    //cp visc for phase 1
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 0.02;
    }

    //  Calculation of properties for phase 2
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real cp_2(const Real Temp1)
    //cp calc for solid pellet phase 2  
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real cp_2a = 0.5867;
      amrex::Real T0 = 373.3;
      amrex::Real T1 = 997.570;
      amrex::Real T2 = 1040.25;
      amrex::Real Tm1= 1007.3125;
      amrex::Real Tm2=Tm1;
      amrex::Real Tm3=1429.125;
      amrex::Real T2A=1016.875;
      amrex::Real T3A=1059.4;
      amrex::Real c1= 2.5489;
      amrex::Real c2= c1;
      amrex::Real c3=8.698;
      amrex::Real o1=9.5625;
      amrex::Real o2=o1;
      amrex::Real o3=10.625;
      amrex::Real cp0=0.5867;
      amrex::Real cp1=cp0;
      amrex::Real cp2=1.1242;
      amrex::Real cp3=0.8144;
      amrex::Real a2 = -0.00729;
      amrex::Real a3 = 0.00014857;
   

      if (Temp<= T0)
	{ cp_2a = cp0;
	}
      else
	{ if (Temp <= T1)
	    { cp_2a = cp1 + c1 * exp(-0.5*pow(((Temp-Tm1)/o1),2));
	    }
	      else
		{ if (Temp <= T1)
		    { cp_2a = cp2 + a2*(Temp-T2A) + c2 * exp(-0.5*pow(((Temp-Tm2)/o2),2));
		    }
		      else
			{  cp_2a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o3),2));
			}}}
   
      return 1000.0*cp_2a;
    }
  
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real k_2(const Real Temp1)
    //k calc for solid pellet phase 2  
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 2.90;
    }
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_2(const Real Temp1)
    //rho calc for solid pellet phase 2  
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 7700.0;
    }

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_2(const Real Temp1)
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_2a = 0.0;
      amrex::Real MeltTemp = 1425.0;
      if (Temp >= MeltTemp)
	{ liqfrac_2a = 1.0;
	}
      else
	{ liqfrac_2a = 0.0;
	}
      return liqfrac_2a;
    }

  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_2(const Real Temp1)
    //cp visc for phase 2
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 0.02;
    }

  
    //  Calculation of properties for phase 3: 1 % carbon bath
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real cp_3(const Real Temp1)
    //cp calc for 1 % C bath
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real cp_3a = 0.6024;
      amrex::Real T0  = 625.0;
      amrex::Real T1  = 1019.0;
      amrex::Real T2  = 1033.875;
      amrex::Real T3  = 1248.5;
      amrex::Real T4  = 1597.0;
      amrex::Real T5  = 1733.0;
      amrex::Real Tm1 = 1026.4375;
      amrex::Real Tm2 = Tm1;
      amrex::Real Tm3 = Tm2;
      amrex::Real Tm4 = 1665.0;
      amrex::Real Tm5 = Tm4;
      amrex::Real Tm6 = Tm5;
      amrex::Real T1A = T0;
      amrex::Real T2A = T1;
      amrex::Real T3A = T2;
      amrex::Real T4A = T3;
      amrex::Real T5A = T4;
      amrex::Real T6A = T5;
      amrex::Real c1  = 2.659163;
      amrex::Real c2  = c1;
      amrex::Real c3  = c1;
      amrex::Real c4  = 2.5421731;
      amrex::Real c5  = c4;
      amrex::Real c6  = c4;
      amrex::Real o1  = 7.4375;
      amrex::Real o2  = o1;
      amrex::Real o3  = o1;
      amrex::Real o4  = 34.0;
      amrex::Real o5 = o4;
      amrex::Real o6 = o5;
      amrex::Real cp0 = .6204;
      amrex::Real cp1 = cp0;
      amrex::Real cp2 = 1.0285;
      amrex::Real cp3 = 0.7207;
      amrex::Real cp4 = 0.6499;
      amrex::Real cp5 = 0.6947;
      amrex::Real cp6 = 0.8683;
      amrex::Real a1 = 0.001081503;
      amrex::Real a2 =-0.0206924;
      amrex::Real a3 =-0.000329848063;
      amrex::Real a4 = 0.00015062494;
      amrex::Real a5 = 0.001276471;
      amrex::Real a6 =-0.0000808824; 
   

      if (Temp<= T0)
	{ cp_3a = cp0;
	}
      else
	{ if (Temp <= T1)
	    { cp_3a = cp1 + a1*(Temp-T1A)+c1 * exp(-0.5*pow(((Temp-Tm1)/o1),2));
	    }
	      else
		{ if (Temp <= T2)
		    { cp_3a = cp2 + a2*(Temp-T2A) + c2 * exp(-0.5*pow(((Temp-Tm2)/o2),2));
		    }
		  else
		    { if (Temp <= T3)
		        { cp_3a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o2),2));
		        }
		      else
		      	{ if (Temp <= T4)
		            { cp_3a = cp4 + a4*(Temp-T4A) + c4 * exp(-0.5*pow(((Temp-Tm4)/o4),2));
		            }
			  else
			    { if (Temp <= T5)
		                { cp_3a = cp5 + a5*(Temp-T5A) + c5 * exp(-0.5*pow(((Temp-Tm5)/o5),2));
		             }
		              else
			        {  cp_3a = cp6 + a6*(Temp-T6A) + c6 * exp(-0.5*pow(((Temp-Tm6)/o6),2));
				}
			        }}}}}
   
      return 1000.0*cp_3a;
    }
  
    
  
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real k_3(const Real Temp1)
    //k calc for 1 % C bath
    {
      amrex::Real Temp = Temp1 + 273.15;
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
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_3(const Real Temp1)
    //rho calc for 1 % C bath
    {
      amrex::Real Temp = Temp1 + 273.15;
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

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_3(const Real Temp1)
    // calculation for 1% carbon bath
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_3a = 0.0;
      amrex::Real MeltTemp = 1669.0;
      if (Temp >= MeltTemp)
	{ liqfrac_3a = 1.0;
	}
      else
	{ liqfrac_3a = 0.0;
	}
      return liqfrac_3a;
    }


  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_3(const Real Temp1)
    //cp visc for phase 3
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 0.02;
    }

  
  
    //  Calculation of properties for phase 4
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real cp_4(const Real Temp1)
    //cp calc for 3% carbon bath  


    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real cp_4a = 0.673;

      amrex::Real T0  = 625.0;
      amrex::Real T1  = 1014.75;
      amrex::Real T2  = 1031.73;
      amrex::Real T3  = 1271.875;
      amrex::Real T4  = 1422.75;
      amrex::Real T5  = 1554.5;
      amrex::Real Tm1 = 1023.25;
      amrex::Real Tm2 = Tm1;
      amrex::Real Tm3 = Tm2;
      amrex::Real Tm4 = 1478.625;
      amrex::Real Tm5 = Tm4;
      amrex::Real Tm6 = Tm5;
      amrex::Real T1A = T0;
      amrex::Real T2A = T1;
      amrex::Real T3A = T2;
      amrex::Real T4A = T3;
      amrex::Real T5A = T4;
      amrex::Real T6A = T5;
      amrex::Real c1  = 2.50463;
      amrex::Real c2  = c1;
      amrex::Real c3  = c1;
      amrex::Real c4  = 2.28955;
      amrex::Real c5  = c4;
      amrex::Real c6  = c4;
      amrex::Real o1  = 8.5;
      amrex::Real o2  = o1;
      amrex::Real o3  = o1;
      amrex::Real o4  = 32.9375;
      amrex::Real o5  = o4;
      amrex::Real o6  = o5;
      amrex::Real cp0 = 0.6073;
      amrex::Real cp1 = cp0;
      amrex::Real cp2 = 0.9792;
      amrex::Real cp3 = 0.8840;
      amrex::Real cp4 = 0.7913;
      amrex::Real cp5 = 0.8599;
      amrex::Real cp6 = 0.8377;
      amrex::Real a1  = 0.0009541416;
      amrex::Real a2  =-0.0056;
      amrex::Real a3  =-0.000386015;
      amrex::Real a4  = 0.0004546050;
      amrex::Real a5  =-0.000168501;
      amrex::Real a6  = 0.0000688525; 
   

      if (Temp<= T0)
	{ cp_4a = cp0;
	}
      else
	{ if (Temp <= T1)
	    { cp_4a = cp1 + a1*(Temp-T1A)+c1 * exp(-0.5*pow(((Temp-Tm1)/o1),2));
	    }
	      else
		{ if (Temp <= T2)
		    { cp_4a = cp2 + a2*(Temp-T2A) + c2 * exp(-0.5*pow(((Temp-Tm2)/o2),2));
		    }
		  else
		    { if (Temp <= T3)
		        { cp_4a = cp3 + a3*(Temp-T3A) + c3 * exp(-0.5*pow(((Temp-Tm3)/o3),2));
		        }
		      else
		      	{ if (Temp <= T4)
		            { cp_4a = cp4 + a4*(Temp-T4A) + c4 * exp(-0.5*pow(((Temp-Tm4)/o4),2));
		            }
			  else
			    { if (Temp <= T5)
		                { cp_4a = cp5 + a5*(Temp-T5A) + c5 * exp(-0.5*pow(((Temp-Tm5)/o5),2));
		             }
		              else
			        {  cp_4a = cp6 + a6*(Temp-T6A) + c6 * exp(-0.5*pow(((Temp-Tm6)/o6),2));
				}
			        }}}}}
   
      return 1000.0*cp_4a;
    }
  
      
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real k_4(const Real Temp1)
    //k calc for 3% carbon bath  
     {
      amrex::Real Temp = Temp1 + 273.15;
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

    //  Calc
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_4(const Real Temp1)
    //rho calc for 3% carbon bath  
     {
      amrex::Real Temp = Temp1 + 273.15;
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
  

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_4(const Real Temp1)
    // calculation for 3% carbon bath
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_4a = 0.0;
      amrex::Real MeltTemp = 1477.0;
      if (Temp >= MeltTemp)
	{ liqfrac_4a = 1.0;
	}
      else
	{ liqfrac_4a = 0.0;
	}
      return liqfrac_4a;
    }


  
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_4(const Real Temp1)
    //cp visc for phase 4
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 0.02;
    }

  

    //  Calculation of properties for phase 5: slag
  
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real cp_5(const Real Temp1)
    //cp calc for slag
       {
      amrex::Real Temp = Temp1 + 273.15;
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
  
    AMREX_GPU_DEVICE AMREX_INLINE
    
    Real k_5(const Real Temp1)
    //k calc for slag  
    {
      amrex::Real Temp = Temp1 + 273.15;
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
   
    AMREX_GPU_DEVICE AMREX_INLINE
   
    Real rho_5(const Real Temp1)
    //rho calc for slag  
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real rho_5a = 2700.0;
      amrex::Real MeltTemp = 1574.0;
      if (Temp >= MeltTemp)
	{ rho_5a = 2700.0;
	}
      else
	{ rho_5a = 3000.0;
	}
      return rho_5a;
    }

    AMREX_GPU_DEVICE AMREX_INLINE

    Real liqfrac_5(const Real Temp1)

    // phase 5 values are for slag  
    // a value of 1 indicates liquid, a value of 0 is solid.
    {
      amrex::Real Temp = Temp1 + 273.15;
      amrex::Real liqfrac_5a = 0.0;
      amrex::Real MeltTemp = 1574.0;
      if (Temp >= MeltTemp)
	{ liqfrac_5a = 1.0;
	}
      else
	{ liqfrac_5a = 0.0;
	}
      return liqfrac_5a;
    }

  
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real visc_5(const Real Temp1)
    //calc visc for phase 5
    {
      amrex::Real Temp = Temp1 + 273.15;
      return 0.02;
    }


  //  Functions to find properties based on temperature and phase
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real get_cp(const Real Temp,int material0)

    // calculate cp based on phase and temp

    {
      Real cpcalc=1000.0;
      int material1 = material0;
      switch (material1) {
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
           cpcalc = cp_3(Temp);
        break;
        case 4:
           cpcalc = cp_4(Temp);
        break;
        case 5:
           cpcalc = cp_5(Temp);
        break;
      }
      return cpcalc;
    }
  
    AMREX_GPU_DEVICE AMREX_INLINE

    Real get_liqfrac(const Real Temp,int material0)

    // calculate liqfrac based on phase and temp

    {
      Real liqfraccalc=0.0;
      int material1 = material0;
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
  
   AMREX_GPU_DEVICE AMREX_INLINE

    Real get_k(const Real Temp,int material0)

    // calculate rho based on phase and temp

    {
      Real kcalc=1000.0;
      int material1 = material0;
      switch (material1) {
	case 0:
           kcalc = k_0(Temp);
	   break;
        case 1:
           kcalc = k_1(Temp);
        break;
        case 2:
           kcalc = k_2(Temp);
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

  
   AMREX_GPU_DEVICE AMREX_INLINE

    Real get_rho(const Real Temp,int material0)

    // calculate rho based on phase and temp

    {
      Real rhocalc=1000.0;
      int material1 = material0;
      switch (material1) {
	case 0:
           rhocalc = rho_0(Temp);
        break;
        case 1:
           rhocalc = rho_1(Temp);
        break;
        case 2:
           rhocalc = rho_2(Temp);
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
     AMREX_GPU_DEVICE AMREX_INLINE

    Real get_visc(const Real Temp,int material0)

    // calculate rho based on phase and temp

    {
      Real visccalc=1000.0;
      int material1 = material0;
      switch (material1) {
	case 0:
           visccalc = visc_0(Temp);
        break;
        case 1:
           visccalc = visc_1(Temp);
        break;
        case 2:
           visccalc = visc_2(Temp);
        break;
        case 3:
           visccalc = visc_3(Temp);
        break;
        case 4:
           visccalc = visc_4(Temp);
        break;
        case 5:
           visccalc = visc_5(Temp);
        break;
      }
      return visccalc;
    }
  

  //   End of properties calculations subroutines 

  
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
	    amrex::Real Tinit_liq,Tinit_solid;

	    int  solidmat,liquidmat;

	    //  Set inital temps and materials, this is temporary code.

        solidmat = 2;
        liquidmat = 5;

        Tinit_liq = 1500.0;

        Tinit_solid = 30.0;

	

        // get solid  properties
        mol_fe = bound01(get_liqfrac(Temp,2));
		
        sol_fe = bound01((1.0 - mol_fe));
				
	
        cp_fe  = get_cp(Temp,solidmat);
        dens_fe = get_rho(Temp,solidmat);
        cond_fe   = get_k(Temp,solidmat);
       
        
        // get slag properties
        mol_slg = bound01(get_liqfrac(Temp,liquidmat));
                    
        sol_slg = bound01(1.0 - mol_slg);

	
	
        cp_slg  = get_cp(Temp,liquidmat);
        dens_slg = get_rho(Temp,liquidmat);
        cond_slg   = get_k(Temp,liquidmat);
            
        amrex::Real vfrac_fe=bound01((phi(i,j,k,NTHERMVARS+MIXMASS_ID)-dens_slg)
				     /(dens_fe-dens_slg));
        phi(i,j,k,DENS_ID)   = dens_slg*(1.0-vfrac_fe)+dens_fe*vfrac_fe;
        phi(i,j,k,COND_ID)   = cond_slg*(1.0-vfrac_fe)+cond_fe*vfrac_fe;
        phi(i,j,k,SPHEAT_ID)   = cp_slg*(1.0-vfrac_fe)+cp_fe*vfrac_fe;
        phi(i,j,k,NTHERMVARS+SOLFE_ID)   =vfrac_fe*sol_fe;             
        phi(i,j,k,NTHERMVARS+MOLFE_ID)   =vfrac_fe*mol_fe;             
        phi(i,j,k,NTHERMVARS+SOLSLG_ID)   = bound01((1.0-vfrac_fe))*sol_slg;             
        phi(i,j,k,NTHERMVARS+MOLSLG_ID)   = bound01((1.0-vfrac_fe))*mol_slg;             
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
