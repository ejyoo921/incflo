#! /bin/bash
set -x  # enables all executed commands are printed to the terminal
paren=`pwd`
incflo="${paren}/incflo3d.gnu.x86-spr.MPI.EB.ex"

res=( 32 ) # Select the res == n_cell you want
for i in "${res[@]}"
do
    rm -rf "${i}"
    mkdir "${i}"
    cd "${i}" || exit
    cp "${paren}/Steel-EB.inp" .

    "${incflo}" Steel-EB.inp > out

    ls -1v *plt*/Header | tee movie.visit
    cd "${paren}" || exit
done
