#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -j y
#$ -N psi4H1
#$ -l h_data=500M,h_rt=96:00:00
#$ -pe mpirun1 4
#$ -R y


LD_LIBRARY_PATH=/opt/acml/acml5.1.0/gfortran64_mp/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

/opt/mackerell/apps/psi/psi -i form_water_mp2_bsse_H1.inp -o form_water_mp2_bsse_H1.out -n $NSLOTS

