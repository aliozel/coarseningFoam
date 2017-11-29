#!/bin/bash

cd $1

for var in "mappedVoidfraction" "mappedCoordNumber" "mappedUs"
do
	cp $var dummy
	nCell=$(sed "1,20d" dummy | head -1 )
	sed "1,22d" dummy | sed -e "s=(==g" -e "s=)==g" | head -"$nCell" > $var""post		
done

paste mappedVoidfractionpost mappedCoordNumberpost | awk '{sum+=(1.-$1)*$2;sum2+=(1.-$1)}END{print sum/sum2,sum2/NR}'
paste mappedVoidfractionpost mappedUspost | awk '{sumX+=(1.-$1)*$2;sumY+=(1.-$1)*$3;sumZ+=(1.-$1)*$4;sum2+=(1.-$1)}END{print sumX/sum2,sumY/sum2,sumZ/sum2,sum2/NR}'

rm mappedVoidfractionpost mappedCoordNumberpost mappedUspost dummy

cd ..
