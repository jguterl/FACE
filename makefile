################################################################################
# Automatically-generated file. Do not edit!
################################################################################

#-include ../makefile.init

RM := rm -rf
USER_OBJS=
LIBS=
# All of the sources participating in the build are defined here
# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
src/face.f90 \
src/modFACE_IO.f90 \
src/modFACE_allocate.f90 \
src/modFACE_cap.f90 \
src/modFACE_compute.f90 \
src/modFACE_coupling.f90 \
src/modFACE_error.f90 \
src/modFACE_functions.f90 \
src/modFACE_header.f90 \
src/modFACE_help.f90 \
src/modFACE_init.f90 \
src/modFACE_input.f90 \
src/modFACE_interface.f90 \
src/modFACE_main.f90 \
src/modFACE_maths.f90 \
src/modFACE_misc.f90 \
src/modFACE_output.f90 \
src/modFACE_parser.f90 \
src/modFACE_precision.f90 \
src/modFACE_run.f90 \
src/modFACE_solverf90.f90 \
src/modFACE_step.f90 

F_SRCS += \
src/modFACE_solver.f 

OBJS += \
./src/face.o \
./src/modFACE_IO.o \
./src/modFACE_allocate.o \
./src/modFACE_cap.o \
./src/modFACE_compute.o \
./src/modFACE_coupling.o \
./src/modFACE_error.o \
./src/modFACE_functions.o \
./src/modFACE_header.o \
./src/modFACE_help.o \
./src/modFACE_init.o \
./src/modFACE_input.o \
./src/modFACE_interface.o \
./src/modFACE_main.o \
./src/modFACE_maths.o \
./src/modFACE_misc.o \
./src/modFACE_output.o \
./src/modFACE_parser.o \
./src/modFACE_precision.o \
./src/modFACE_run.o \
./src/modFACE_solver.o \
./src/modFACE_solverf90.o \
./src/modFACE_step.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: src/%.f90
	@echo 'Building file: $<'
	gfortran -funderscoring -O3  -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/face.o: src/face.f90 src/modFACE_interface.o src/modFACE_main.o

src/modFACE_IO.o: src/modFACE_IO.f90 src/modFACE_allocate.o src/modFACE_error.o src/modFACE_header.o src/modFACE_step.o

src/modFACE_allocate.o: src/modFACE_allocate.f90 src/modFACE_header.o src/modFACE_output.o

src/modFACE_cap.o: src/modFACE_cap.f90 src/modFACE_header.o

src/modFACE_compute.o: src/modFACE_compute.f90 src/modFACE_cap.o src/modFACE_error.o src/modFACE_functions.o src/modFACE_header.o

src/modFACE_coupling.o: src/modFACE_coupling.f90 src/modFACE_error.o src/modFACE_header.o

src/modFACE_error.o: src/modFACE_error.f90 src/modFACE_header.o src/modFACE_precision.o

src/modFACE_functions.o: src/modFACE_functions.f90 src/modFACE_error.o src/modFACE_header.o

src/modFACE_header.o: src/modFACE_header.f90 src/modFACE_precision.o

src/modFACE_help.o: src/modFACE_help.f90 src/modFACE_IO.o src/modFACE_misc.o

src/modFACE_init.o: src/modFACE_init.f90 src/modFACE_IO.o src/modFACE_error.o src/modFACE_functions.o src/modFACE_header.o src/modFACE_input.o src/modFACE_output.o src/modFACE_precision.o src/modFACE_step.o

src/modFACE_input.o: src/modFACE_input.f90 src/modFACE_allocate.o src/modFACE_header.o src/modFACE_output.o src/modFACE_parser.o

src/modFACE_interface.o: src/modFACE_interface.f90 src/modFACE_coupling.o src/modFACE_error.o src/modFACE_header.o src/modFACE_main.o src/modFACE_precision.o

src/modFACE_main.o: src/modFACE_main.f90 src/modFACE_IO.o src/modFACE_header.o src/modFACE_help.o src/modFACE_run.o

src/modFACE_maths.o: src/modFACE_maths.f90 src/modFACE_error.o src/modFACE_header.o src/modFACE_precision.o

src/modFACE_misc.o: src/modFACE_misc.f90 src/modFACE_output.o

src/modFACE_output.o: src/modFACE_output.f90 src/modFACE_error.o src/modFACE_header.o

src/modFACE_parser.o: src/modFACE_parser.f90 src/modFACE_header.o src/modFACE_help.o src/modFACE_misc.o src/modFACE_output.o

src/modFACE_precision.o: src/modFACE_precision.f90

src/modFACE_run.o: src/modFACE_run.f90 src/modFACE_coupling.o src/modFACE_header.o src/modFACE_init.o src/modFACE_input.o

src/%.o: src/%.f
	@echo 'Building file: $<'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/modFACE_solver.o: src/modFACE_solver.f src/modFACE_compute.o src/modFACE_error.o src/modFACE_functions.o src/modFACE_header.o src/modFACE_output.o src/modFACE_precision.o src/modFACE_solverf90.o

src/modFACE_solverf90.o: src/modFACE_solverf90.f90 src/modFACE_compute.o src/modFACE_error.o src/modFACE_functions.o src/modFACE_header.o src/modFACE_maths.o src/modFACE_output.o src/modFACE_precision.o

src/modFACE_step.o: src/modFACE_step.f90 src/modFACE_compute.o src/modFACE_functions.o src/modFACE_header.o src/modFACE_output.o src/modFACE_solver.o


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: FACE2.0

# Tool invocations
FACE2.0: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	gfortran  -o "FACE2.0" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) FACE2.0
	-@echo ' '

.PHONY: all clean dependents

