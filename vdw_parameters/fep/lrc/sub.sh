#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -o dynamics.txt
#$ -j y
#$ -N LRC
#$ -l h_data=200M,h_rt=196:00:00
#$ -pe mpirun1  8 
#$ -R y
#  -q all.q@compute-0-76

charmm=/opt/mackerell/apps/charmm/serial/c38b2-serial

#Path to most appropriate mpirun is set system-wide by Ron,
#but we still need --no-daemonize
mpirun="mpirun --leave-session-attached"
ulimit -c 4

$mpirun -np $NSLOTS $charmm -i run_box_md.inp > run_box_md.out

$charmm < anal_box_md_lrc.inp > anal_box_md_lrc.out 

