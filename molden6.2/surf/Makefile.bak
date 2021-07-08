# Compilation flags
#CC	    = cc
CC	    = cc
INCLUDE     = -I. 
#LINCLUDE    = -lcurses -ltermcap -lm
LINCLUDE    = -lm
OPT_CFLAGS  = -O2 $(FLAGS) $(INCLUDE)
#CFLAGS	    = -g $(FLAGS) $(INCLUDE)
CFLAGS	    = -O2 $(FLAGS) $(INCLUDE) ${EXT} -Wno-implicit-function-declaration

ifeq ($(uname), Darwin)
#
# Mac OS X g77
#
CFLAGS= -std=gnu89 -Wno-return-type -DDARWIN -O2 $(FLAGS) $(INCLUDE) -Wno-implicit-function-declaration
endif

# These are the user object files in the application
SRCS        = surf.c io.c compute.c dual.c utils.c lp.c chull.c tessel_cases.c \
	      tessel_patches.c tessel_convex.c tessel_concave.c tessel_torus.c 

# obj files for sequential execution
OBJS        = surf.o io.o compute.o dual.o utils.o lp.o chull.o tessel_cases.o \
	      tessel_patches.o tessel_convex.o tessel_concave.o tessel_torus.o  

# make objects
surf: $(OBJS) Makefile 
	$(CC) $(CFLAGS) $(OBJS) -o surf $(LINCLUDE)

lint:
	lint $(INCLUDE) $(SRCS)

clean:
	-rm -f *.o 

release:
	tar -cvf surf.tar README *.[hc] Makefile
	compress surf.tar

.c.o: 
	$(CC) $(CFLAGS) -c $*.c


# make depend makes the proper include file dependencies.  You _could_ run
# it on a sun4, but there's a bug in the SunOS version of sed that causes
# the all-important "cat makedep | sed -e 's/\.o/.d4o/p'" etc. line to
# not print the results to the stdout (but to keep pumping it through the
# edits.  The 'p' flag effectively does nothing). To avoid this I have
# replaced the 'p' flag by the 'w file' flag in sed.

DEPEND = $(SRCS)

depend: $(DEPEND)
	@ echo making dependencies...
	@ echo ' ' > makedep
	@ makedepend $(INCLUDE) -f makedep $(DEPEND)
	@ echo '/^# DO NOT DELETE THIS LINE/+2,$$d' >eddep
	@ echo '$$r makedep' >>eddep
	@ echo 'w' >>eddep
	@ cp Makefile Makefile.bak
	@ ex - Makefile < eddep
	@ rm eddep makedep 
	@ echo '  ' >> Makefile
	@ echo '# DEPENDENCIES MUST END AT END OF FILE' >> Makefile
	@ echo '# IF YOU PUT STUFF HERE IT WILL GO AWAY' >> Makefile
	@ echo '# see make depend above' >> Makefile
	@ echo 'done.'
	@ echo '  ' >> Makefile

# DO NOT DELETE THIS LINE -- make depend depends on it.
  
 
# DO NOT DELETE

surf.o: surf.h linalg.h
io.o: surf.h linalg.h
compute.o: surf.h linalg.h chull.h dual.h
dual.o: surf.h linalg.h dual.h chull.h
utils.o: surf.h linalg.h
lp.o: surf.h linalg.h
chull.o: surf.h linalg.h chull.h
tessel_cases.o: surf.h linalg.h dual.h
tessel_patches.o: surf.h linalg.h
tessel_convex.o: surf.h linalg.h
tessel_concave.o: surf.h linalg.h
tessel_torus.o: surf.h linalg.h
  
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
  
