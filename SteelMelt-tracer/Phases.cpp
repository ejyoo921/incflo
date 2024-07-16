#include<Phases.H>

namespace phasephysics
{
    amrex::Vector<std::string> phasenames(NUM_PHASES);

    void init()
    {
        phasenames[SOLFE_ID]="Solid_Fe";
        phasenames[MOLFE_ID]="Molten_Fe";
        phasenames[SOLSLG_ID]="Solid_slag";
        phasenames[MOLSLG_ID]="Molten_slag";
        phasenames[MIXMASS_ID]="Mixture_mass";
        phasenames[VFRAC_ID]="volume_fraction_fe";
    }    
    
    void close()
    {
        phasenames.clear();
    }
    
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(phasenames.begin(),phasenames.end(),specname);
        if(it != phasenames.end())
        {
            loc=it-phasenames.begin();
        }
        return(loc);
    }
}
