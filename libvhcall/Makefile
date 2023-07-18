#FC = /opt/nec/ve/bin/nfort
FC = mpinfort
GFORT = gfortran
#INC_VH=-I ../build_vh
#INC_VE=-I ../build_ve
ALL = libvhcall.so #sample
all: $(ALL)
sample: sample.f90
	$(FC) $(INC_VE) -o $@ $^ -lvhcall_fortran
libvhcall.so: libvhcall.f90
	$(GFORT) $(INC_VH) -o $@ $^ -fPIC -shared
clean:
	rm -f $(ALL) *.o *.mod
