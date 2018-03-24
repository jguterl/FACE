################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../face.f90 \
../modFACE_IO.f90 \
../modFACE_allocate.f90 \
../modFACE_cap.f90 \
../modFACE_compute.f90 \
../modFACE_coupling.f90 \
../modFACE_error.f90 \
../modFACE_functions.f90 \
../modFACE_header.f90 \
../modFACE_help.f90 \
../modFACE_init.f90 \
../modFACE_input.f90 \
../modFACE_interface.f90 \
../modFACE_main.f90 \
../modFACE_maths.f90 \
../modFACE_misc.f90 \
../modFACE_output.f90 \
../modFACE_parser.f90 \
../modFACE_precision.f90 \
../modFACE_run.f90 \
../modFACE_solverf90.f90 \
../modFACE_step.f90 

F_SRCS += \
../modFACE_solver.f 

OBJS += \
./face.o \
./modFACE_IO.o \
./modFACE_allocate.o \
./modFACE_cap.o \
./modFACE_compute.o \
./modFACE_coupling.o \
./modFACE_error.o \
./modFACE_functions.o \
./modFACE_header.o \
./modFACE_help.o \
./modFACE_init.o \
./modFACE_input.o \
./modFACE_interface.o \
./modFACE_main.o \
./modFACE_maths.o \
./modFACE_misc.o \
./modFACE_output.o \
./modFACE_parser.o \
./modFACE_precision.o \
./modFACE_run.o \
./modFACE_solver.o \
./modFACE_solverf90.o \
./modFACE_step.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	/usr/local/bin/gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

face.o: ../face.f90 modFACE_interface.o modFACE_main.o

modFACE_IO.o: ../modFACE_IO.f90 modFACE_allocate.o modFACE_error.o modFACE_header.o modFACE_step.o

modFACE_allocate.o: ../modFACE_allocate.f90 modFACE_header.o modFACE_output.o

modFACE_cap.o: ../modFACE_cap.f90 modFACE_header.o

modFACE_compute.o: ../modFACE_compute.f90 modFACE_cap.o modFACE_error.o modFACE_functions.o modFACE_header.o

modFACE_coupling.o: ../modFACE_coupling.f90 modFACE_error.o modFACE_header.o

modFACE_error.o: ../modFACE_error.f90 modFACE_header.o modFACE_precision.o

modFACE_functions.o: ../modFACE_functions.f90 modFACE_error.o modFACE_header.o

modFACE_header.o: ../modFACE_header.f90 modFACE_precision.o

modFACE_help.o: ../modFACE_help.f90 modFACE_IO.o modFACE_misc.o

modFACE_init.o: ../modFACE_init.f90 modFACE_IO.o modFACE_error.o modFACE_functions.o modFACE_header.o modFACE_input.o modFACE_output.o modFACE_precision.o modFACE_step.o

modFACE_input.o: ../modFACE_input.f90 modFACE_allocate.o modFACE_header.o modFACE_output.o modFACE_parser.o

modFACE_interface.o: ../modFACE_interface.f90 modFACE_coupling.o modFACE_error.o modFACE_header.o modFACE_main.o modFACE_precision.o

modFACE_main.o: ../modFACE_main.f90 modFACE_IO.o modFACE_header.o modFACE_help.o modFACE_run.o

modFACE_maths.o: ../modFACE_maths.f90 modFACE_error.o modFACE_header.o modFACE_precision.o

modFACE_misc.o: ../modFACE_misc.f90 modFACE_output.o

modFACE_output.o: ../modFACE_output.f90 modFACE_error.o modFACE_header.o

modFACE_parser.o: ../modFACE_parser.f90 modFACE_header.o modFACE_help.o modFACE_misc.o modFACE_output.o

modFACE_precision.o: ../modFACE_precision.f90

modFACE_run.o: ../modFACE_run.f90 modFACE_coupling.o modFACE_header.o modFACE_init.o modFACE_input.o

%.o: ../%.f
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	/usr/local/bin/gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

modFACE_solver.o: ../modFACE_solver.f modFACE_compute.o modFACE_error.o modFACE_functions.o modFACE_header.o modFACE_output.o modFACE_precision.o modFACE_solverf90.o

modFACE_solverf90.o: ../modFACE_solverf90.f90 modFACE_compute.o modFACE_error.o modFACE_functions.o modFACE_header.o modFACE_maths.o modFACE_output.o modFACE_precision.o

modFACE_step.o: ../modFACE_step.f90 modFACE_compute.o modFACE_functions.o modFACE_header.o modFACE_output.o modFACE_solver.o


