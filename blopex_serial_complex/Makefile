objects = multivector pcg_multi matmultivec

driver: $(objects) FORCE
	cd driver && $(MAKE)	
        
multivector: FORCE
	cd multivector && $(MAKE)
        
pcg_multi: FORCE
	cd pcg_multi && $(MAKE)
        
matmultivec: FORCE
	cd matmultivec && $(MAKE)

clean: FORCE
	cd multivector && $(MAKE) clean
	cd pcg_multi && $(MAKE) clean
	cd driver && $(MAKE) clean
	cd matmultivec && $(MAKE) clean

FORCE:
