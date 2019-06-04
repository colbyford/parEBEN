#!/bin/sh

#PBS -q steelhead
#PBS -N parEBEN_doMPI_FullTest
#PBS -l nodes=6:ppn=36,mem=2000GB
#PBS -l walltime=200:00:00

module load openmpi/1.10.0 R/3.5.0
cd $PBS_O_WORKDIR

nodecnt=`wc -l $PBS_NODEFILE | awk '{print $1}'`
echo "Node File:"
cat $PBS_NODEFILE
echo ""

echo "Executing: mpirun -quiet -np 1 --hostfile $PBS_NODEFILE R --slave -f Full_Test_19872_Gaus_doMPI.R"
echo ""
mpirun -quiet -np 1 --hostfile $PBS_NODEFILE R --slave -f Full_Test_19872_Gaus_doMPI.R
echo ""
