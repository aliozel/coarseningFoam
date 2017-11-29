#!/bin/bash

awk '{sumalppV+=$15*(1.-$18);sumalpp+=(1.-$18)}END{print "alpp =",sumalpp/NR, "Us =",sumalppV/sumalpp}' dummy

mp=$(awk 'BEGIN{print 1./6.*3.14152*75e-06*75e-06*75e-06*1500}')
taupSt="0.025"
Volmesh="1.3824e-11"

awk -v mp=$mp -v taupSt=$taupSt -v Volmesh=$Volmesh '{sumdrag+=mp/taupSt*(0.0696544-$15)*$21}END{print "Fdrag =",sumdrag/Volmesh}' dummy
#awk -v mp=$mp -v taupSt=$taupSt -v Volmesh=$Volmesh '{print mp/taupSt*($10-$15)/Volmesh*$21}' dummy
