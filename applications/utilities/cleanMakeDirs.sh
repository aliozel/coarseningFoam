#/bin/bash

for folder in */; 
do 
	cd ./$folder
	wclean
	rm -rf .git
	if [ -d Make ]; then
	   cd ./Make 		 
	   echo $PWD
	   rm -rf linux*
	   rm -rf darwin*	
	   cd ..
	fi 
	cd ..
done
