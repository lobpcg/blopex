driver: petsc_interface
	cd driver && $(MAKE)
	cd driver_fiedler && $(MAKE)
        
petsc_interface: FORCE
	cd petsc-interface && $(MAKE)

clean: FORCE
	cd petsc-interface && $(MAKE) clean
	cd driver && $(MAKE) clean
	cd driver_fiedler && $(MAKE) clean

FORCE:
