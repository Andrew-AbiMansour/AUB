VERSION = 3.0.2

C:
	( cd Lib ; $(MAKE) )
	( cd cv ; $(MAKE) )

all: C cov

library:
	( cd Lib ; $(MAKE) )
	( cd cv; $(MAKE) )
clean:
	( cd Lib ; $(MAKE) clean )

purge:
	( cd Lib ; $(MAKE) purge )

distclean: purge

# do not install CSparse; use CXSparse instead
install:

# uninstall CSparse: do nothing
uninstall:

