PROG = gems.exe

OBJS =	gems_bound.o gems_constant.o gems_data.o gems_disc.o gems_fv.o gems_geom.o \
	gems_input.o gems_jacob.o gems_library.o gems_linear.o gems_main.o \
	gems_output.o gems_precon.o gems_state.o gems_type.o gems_turb.o gems_reaction.o \
	gems_source.o gems_rptbdb.o

#CC = cc
#CFLAGS = -O
#FC = lf95
#FFLAGS = -O --trap --tpp

F90 = mpif90
F90FLAGS = -O0 -g -trace -fpe0
#LDFLAGS = --static

all: $(PROG)

all_clean: $(PROG) clean

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) *.mod $(PROG)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

gems_source.o: gems_state.o
gems_data.o: gems_type.o gems_rptbdb.o
gems_disc.o: gems_state.o
gems_turb.o: gems_state.o
gems_fv.o: gems_disc.o gems_jacob.o gems_precon.o gems_turb.o
gems_bound.o: gems_disc.o
gems_geom.o: gems_data.o gems_library.o
gems_input.o: gems_bound.o gems_geom.o gems_rptbdb.o gems_turb.o gems_reaction.o
gems_reaction.o: gems_jacob.o
gems_jacob.o: gems_state.o
gems_library.o: gems_type.o
gems_linear.o: gems_bound.o
gems_main.o: gems_input.o gems_linear.o gems_output.o gems_source.o
gems_output.o: gems_fv.o gems_reaction.o
gems_precon.o: gems_state.o
gems_state.o: gems_data.o gems_library.o gems_rptbdb.o
gems_type.o: gems_constant.o





