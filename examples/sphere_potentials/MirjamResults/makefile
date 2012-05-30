FCOMPFLAGS    =	-O -ffixed-line-length-132
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	g95
FC            =	g95
CC            =	cc

MAKEFILE      =	makefile
PROGRAM       =	DNAPairPotentials

INCS	      =	

OBJS	      =	DNAPairPotentials.o
		
	
all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
