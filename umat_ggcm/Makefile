clean:
	@rm *.so *.o 2>/dev/null
compile:
	gfortran -c -fPIC -g -O0 TRIPUMAT-HT_5_8_GrainCluster.f
	f2py fortPyBind.f90 -m pyUMATLink -h fortPyBind.pyf --latex-doc --verbose
	f2py -c fortPyBind.pyf fortPyBind.f90 TRIPUMAT-HT_5_8_GrainCluster.o

