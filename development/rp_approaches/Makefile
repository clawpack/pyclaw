PYCLAWMAKE = $(PYCLAW)/Makefile.common

#RP_SOURCE =   rpn2_euler_5wave.f rpt2_euler_5wave.f
#RP_SOURCE =   rpn2_euler_5wave_recompute.f rpt2_euler_5wave_recompute.f
RP_SOURCE =   rpn2_euler_5wave_rec_loc.f rpt2_euler_5wave_rec_loc.f
#RP_SOURCE =   rpn2_euler_5wave_aux.f rpt2_euler_5wave_aux.f
#RP_SOURCE =   rp2_euler_5wave_module.f90

#override TWO_D_CLASSIC_SOURCES = $(TWO_D_CLASSIC)/step2.f $(TWO_D_CLASSIC)/step2ds.f ./flux2.f90 $(TWO_D_CLASSIC)/limiter.f $(TWO_D_CLASSIC)/philim.f

all: classic2.so sharpclaw2.so

include $(PYCLAWMAKE)

