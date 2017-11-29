#!/bin/bash

cd $1

for var in "mappedVoidfraction" "mappedCoordNumber"
do
	cp $var dummy
	nCell=$(sed "1,20d" dummy | head -1 )
	sed "1,22d" dummy | sed -e "s=(==g" -e "s=)==g" | head -"$nCell" > $var""post		
done

paste mappedVoidfractionpost mappedCoordNumberpost | awk '{sum+=(1.-$1)*$2;sum2+=(1.-$1)}END{print sum/sum2}'
rm mappedVoidfractionpost mappedCoordNumberpost dummy

cd ..
