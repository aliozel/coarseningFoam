#!/bin/bash

dp=250e-06
Lx=0.003
Ly=0.003
Lz=0.096

nx=2
ny=$(awk -v nx=$nx -v Lx=$Lx -v Ly=$Ly 'BEGIN{print int(nx*Ly/Lx)}')
nz=$(awk -v nx=$nx -v Lx=$Lx -v Lz=$Lz 'BEGIN{print int(nx*Lz/Lx)}')

echo "Number of particles = [$nx":"$ny":"$nz]"

deltax=$(awk -v Lx=$Lx -v nx=$nx 'BEGIN{print Lx/nx}')
deltay=$(awk -v Ly=$Ly -v ny=$ny 'BEGIN{print Ly/ny}')
deltaz=$(awk -v Lz=$Lz -v nz=$nz 'BEGIN{print Lz/nz}')

echo "Distance = $deltax,$deltay,$deltaz"

nP=$(awk -v nx=$nx -v ny=$ny -v nz=$nz 'BEGIN{print nx*ny*nz}')	
echo $nP > latticePart.dat
for k in $(seq 1 $nz)
do
	for j in $(seq 1 $ny)
	do
		for i in $(seq 1 $nx)
		do			
			xi=$(awk -v i=$i -v deltax=$deltax 'BEGIN{print (i-0.5)*deltax}')
			yi=$(awk -v j=$j -v deltay=$deltay 'BEGIN{print (j-0.5)*deltay}')
			zi=$(awk -v k=$k -v deltaz=$deltaz 'BEGIN{print (k-0.5)*deltaz}')
			echo $xi $yi $zi >> latticePart.dat
		done
	done	
done
