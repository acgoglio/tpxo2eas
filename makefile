ARCH = $(shell uname -s)
ifeq ($(ARCH),Linux)
 FC = ifort
 NCLIB = /users/home/opt-intel_2015.3.187/netcdf/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1/lib/
 NCINCLUDE = /users/home/opt-intel_2015.3.187/netcdf/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1/include/
 NCLIBS= -lnetcdf -lnetcdff
endif

predict_tide: predict_tide.f90 subs.f90 constit.h
	$(FC) -o predict_tide predict_tide.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE) 
	#rm *.o
extract_HC:  extract_HC.f90 subs.f90
	$(FC) -o  extract_HC extract_HC.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE)
	#rm *.o
