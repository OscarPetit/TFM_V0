# Makefile — stable build flags
# Use a generic architecture to avoid illegal instructions
FC      = gfortran

# Remove aggressive flags that may produce unsupported instructions
# Use moderate optimization and debugging info for stability
FFLAGS  = -O3 \
          -g \
          -fbacktrace \
          -fcheck=all \

# List all your source files
SRCS    = 2d_main_final1.f90 boundxy.f90 calculo_funciones_flujos.f90 \
          derivadasx.f90 derivadasy.f90 flujo_yR.f90 get_primitives.f90 \
          initial_conditions.f90 flujo_xR.f90 weno_2d_Lx.f90\
          get_psidot.f90 get_uxdot.f90 get_uydot.f90 \
		  weno_2d_Ly.f90 weno_2d_Rx.f90 weno_2d_Ry.f90

# Generate object file names
OBJS    = $(SRCS:.f90=.o)

# Default target: build the executable "main"
all: main

main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Compile each Fortran source into object
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean up
clean:
	rm -f $(OBJS) main
