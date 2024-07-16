#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

#define SOLFE_ID 0
#define MOLFE_ID 1
#define SOLSLG_ID 2
#define MOLSLG_ID 3
#define MIXMASS_ID 4
#define VFRAC_ID 5
#define NUM_PHASES 6

using namespace amrex;

namespace phasephysics
{
    extern amrex::Vector<std::string> phasenames;
    void init();
    void close();
    int find_id(std::string specname);
}
#endif
