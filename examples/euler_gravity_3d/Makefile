
F2PY = f2py
RIEMANN_SOURCE = ../../../riemann/src/

all:
	$(MAKE) euler_3d_gmap.so
	$(MAKE) mappedGrid.so

euler_3d_gmap.so: $(RIEMANN_SOURCE)/rpn3_euler_mapgrid.f90 $(RIEMANN_SOURCE)/rpt3_euler_mapgrid.f90 $(RIEMANN_SOURCE)/rptt3_euler_mapgrid.f90
	f2py -m $(basename $(notdir $@)) -c $^

mappedGrid.so: $(RIEMANN_SOURCE)/euler3d_mappedGrid.f90
	f2py -m $(basename $(notdir $@)) -c $^

clean:
	rm -f *.so *~ *.pyc
