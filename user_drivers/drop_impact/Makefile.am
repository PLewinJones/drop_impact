# Name of executable
noinst_PROGRAMS = drop_impact

drop_impact_SOURCES = drop_impact.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
drop_impact_LDADD = -L@libdir@ -ldrop_impact_elements -lnavier_stokes -lfluid_interface \
	                     -laxisym_navier_stokes -lsolid -lconstitutive -lgeneric \
			     $(EXTERNAL_LIBS) $(FLIBS)
			     
AM_CPPFLAGS += -I@includedir@
