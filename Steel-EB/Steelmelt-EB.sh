#! /bin/bash
set -x  # enables all executed commands are printed to the terminal
paren=`pwd`
pelec="${paren}/incflo3d.gnu.MPI.EB.ex"
mpi_ranks=8

# pelec_nompi="${paren}/incflo3d.gnu.EB.ex"

res=( 32 ) # Select the res == n_cell you want
for i in "${res[@]}"
do
    rm -rf "${i}"
    mkdir "${i}"
    cd "${i}" || exit
    cp "${paren}/Steel-EB.inp" .

    # Select either one, depending on MPI usage
    mpirun -n ${mpi_ranks} "${pelec}" Steel-EB.inp amr.n_cell="${i} ${i} ${i}" > out
    # "${pelec_nompi}" Steel-EB.inp amr.n_cell="${i} ${i} ${i}" > out

    ls -1v *plt*/Header | tee movie.visit
    cd "${paren}" || exit
done
