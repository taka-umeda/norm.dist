FC              = mpiifort -mcmodel=large -ipo -ip
OPT             = -O3 -xICELAKE-SERVER
LINKER          = $(FC) $(OPT)

PROGRAM         = a.out
MODU		= functions.o
MAIN		= main.o

all:            $(PROGRAM)

$(MODU):	%.o : %.f90
		$(FC) $(OPT) -c $<

$(MAIN):	%.o : %.f90 $(MODU)
		$(FC) $(OPT) -c $<

$(PROGRAM): 	$(MAIN) $(MODU)
		$(LINKER) -o $@ $(MAIN) $(MODU)
