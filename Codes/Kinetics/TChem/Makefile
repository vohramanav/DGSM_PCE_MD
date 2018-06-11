SHELL = /bin/sh

include MakeIncl.Vars

all: library $(EXAMPLEDIR)

library:
	cd src; make; cd ..

ex: 
	cd example/ign-c;    make; cd ../..
clean:
	cd src;     make clean; cd ..
	cd example; make clean; cd ..

distclean: clean
	cd lib;     make distclean; cd ..
	cd example; make distclean; cd .. 


