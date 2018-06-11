include ../MakeIncl.Vars

CVODEARCH = cvode-2.7.0
CVODEINST = $(CVODEARCH)-install
DVODEARCH = dvode
DVODELIB  = lib$(DVODEARCH).a

AR = ar cr

default: cv dv

cv:
	if test ! -d build ; then \
	  mkdir build; \
	else \
	  cd build; /bin/rm -rf *; cd ..;\
	fi
	if test ! -d $(CVODEINST) ; then \
          gunzip $(CVODEARCH).tar.gz; tar xvf $(CVODEARCH).tar ; \
          cd build; cmake ../$(CVODEARCH) -DCMAKE_INSTALL_PREFIX=../$(CVODEARCH)-install; \
	  make ; make install ; cd .. ; gzip $(CVODEARCH).tar;\
	fi

dv:
	if test ! -e $(DVODELIB) ; then \
          gunzip $(DVODEARCH).tar.gz; tar xvf $(DVODEARCH).tar ; \
	  cd $(DVODEARCH); $(F77) $(FFLAGS) -c dvode.f; \
	  $(AR) ../$(DVODELIB) dvode.o; \
	  cd .. ; gzip $(DVODEARCH).tar;\
	fi

clean : distclean

distclean :
	cd build; /bin/rm -rf *; cd ..
	/bin/rm -rf $(CVODEARCH) $(CVODEINST)
	/bin/rm -rf $(DVODEARCH) $(DVODELIB)
	/bin/rm libtchem.a libtchemutil.a

