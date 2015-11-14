#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:19:00
#PBS -j oe
 
cd /gpfs/home/brt10/www/bam2_4 
./bam.x  ./test1.input    
