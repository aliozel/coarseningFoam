//- Filter sizes	  
if(fluidCoarsening) 
{
 //- Fluid Coarsening
 forAll(filterWidth,fWidth)
 {
   //- Create name for filter width			  
   char charName[100];
   sprintf(charName, "%dX%dX%d",filterWidth[fWidth],filterWidth[fWidth],filterWidth[fWidth]);
   word filterWidthName(charName);	

   // Create header for output files
   OFstream* outputFile;
   fileName outpuFilename(mesh.time().path()/outputRelativePath/"pressureGradientBudgetAnalysisZ"+filterWidthName+".dat");
   outputFile =  new OFstream(outpuFilename);
   *outputFile  << "#Time " << "NFirstTerm "  << tab << "RStar" << tab << "N" << tab << "N1"<< "N2"<< tab << "N3" << tab << "barR1"<< tab<< "barR2"<< tab<< "barR3"<< tab<< "barR4"<< tab<< "barR5"<< tab<< "RStarres"<< tab << "barRhoMixVPrime"<< tab
			    << endl;  
   //-Time loop  
   forAll(timeDirs, timeI)
   {
	// Write results into a file
	*outputFile	<< timeDirs[timeI] 							<< tab 
			<< domainAveVars[timeI][fWidth][0] 					<< tab  
                	<< domainAveVars[timeI][fWidth][1]                 			<< tab
	  		<< domainAveVars[timeI][fWidth][2]					<< tab
	  		<< domainAveVars[timeI][fWidth][3]					<< tab
	  		<< domainAveVars[timeI][fWidth][4]					<< tab
	  		<< domainAveVars[timeI][fWidth][5]					<< tab
	  		<< domainAveVars[timeI][fWidth][6]					<< tab
	  		<< domainAveVars[timeI][fWidth][7]					<< tab
	  		<< domainAveVars[timeI][fWidth][8]					<< tab
	  		<< domainAveVars[timeI][fWidth][9]					<< tab
	  		<< domainAveVars[timeI][fWidth][10]					<< tab
	  		<< domainAveVars[timeI][fWidth][11]					<< tab
	  		<< domainAveVars[timeI][fWidth][12]					<< tab


 	<< endl;
   }
  }
 }   
