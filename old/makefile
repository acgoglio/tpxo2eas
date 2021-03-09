ARCH = $(shell uname -s)
ifeq ($(ARCH),Linux)
 FC = ifort
 NCLIB = /zeus/opt/intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/lib/
 NCINCLUDE = /zeus/opt/intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/include/
 NCLIBS= -lnetcdf -lnetcdff
endif

predict_tide: predict_tide.f90 subs.f90 constit.h
	$(FC) -o predict_tide predict_tide.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE) 
	#rm *.o
extract_HC:  extract_HC.f90 subs.f90
	$(FC) -o  extract_HC extract_HC.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE)
	#rm *.o
