int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				

Info<< " " << endl;
particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
Info<< tab <<"Reading particle velocities" << endl;
Info<< " " << endl;

particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
Info<< tab <<"Reading particle positions" << endl;
Info<< " " << endl;

particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
Info<< tab <<"Reading particle radius" << endl;		
Info<< " " << endl;

particleCloud.dataExchangeM().getData("type","scalar-atom",types,count);
Info<< "        Reading particle types " << endl;
Info<< " " << endl;

particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
particleCloud.setPos(positions);
particleCloud.setVel(velocities);
particleCloud.setType(types);		

if(verbose)
{
	forAll(exList,ii)
	{
		int index = exList[ii];
		Pout << tab << "index  = " << index << endl;
		Pout << tab << "rp     = " << particleCloud.radius(index) << endl;
		Pout << tab << "Type   = " << particleCloud.type(index) << endl;
		Pout << tab << "Vp     = " << particleCloud.velocity(index) << endl;
		Pout << tab << "Xp     = " << particleCloud.position(index) << endl;
		Pout << tab << "CellID = " << particleCloud.particleCell(index) << endl;
		if( particleCloud.particleCell(index) > -1 ) Pout << tab << "Uf     = " << U[particleCloud.particleCell(index)] << endl;
		Pout << " " << endl;
	}	
}


			
