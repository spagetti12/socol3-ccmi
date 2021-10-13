LIB =	lapack

SRCS =  $(shell ls *.f)

OBJS =  $(SRCS:.f=.o)

all: $(LIB)

$(LIB): ../lib/lib$(LIB).a 

../lib/lib$(LIB).a: $(OBJS)
	$(AR) $(ARFLAGS) ../lib/lib$(LIB).a $(OBJS) 

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS)



