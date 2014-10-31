euler_3d_gmap.so: rpn3_euler_mapgrid.f90 rpt3_euler_mapgrid.f90 rptt3_euler_mapgrid.f90
	f2py -m $(basename $(notdir $@)) -c $^

mappedGrid.so: euler3d_mappedGrid.f90
	f2py -m $(basename $(notdir $@)) -c $^

clean:
	rm -f *.so *~ *.pyc
