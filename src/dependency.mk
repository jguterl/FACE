$(OBJDIR)/face.o:  $(OBJDIR)/modFACE_interface.o $(OBJDIR)/modFACE_main.o

$(OBJDIR)/modFACE_IO.o:  $(OBJDIR)/modFACE_allocate.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_allocate.o:  $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_output.o

$(OBJDIR)/modFACE_cap.o:  $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_compute.o:  $(OBJDIR)/modFACE_cap.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_functions.o $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_coupling.o: $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_error.o:  $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_precision.o

$(OBJDIR)/modFACE_functions.o:  $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_header.o:  $(OBJDIR)/modFACE_precision.o

$(OBJDIR)/modFACE_help.o:  $(OBJDIR)/modFACE_IO.o $(OBJDIR)/modFACE_misc.o

$(OBJDIR)/modFACE_init.o:  $(OBJDIR)/modFACE_IO.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_functions.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_input.o $(OBJDIR)/modFACE_output.o $(OBJDIR)/modFACE_precision.o $(OBJDIR)/modFACE_step.o

$(OBJDIR)/modFACE_input.o:  $(OBJDIR)/modFACE_allocate.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_output.o $(OBJDIR)/modFACE_parser.o

$(OBJDIR)/modFACE_interface.o:  $(OBJDIR)/modFACE_coupling.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_main.o $(OBJDIR)/modFACE_precision.o

$(OBJDIR)/modFACE_main.o:  $(OBJDIR)/modFACE_IO.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_help.o $(OBJDIR)/modFACE_run.o

$(OBJDIR)/modFACE_maths.o:  $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_precision.o

$(OBJDIR)/modFACE_misc.o:  $(OBJDIR)/modFACE_output.o

$(OBJDIR)/modFACE_output.o: $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_header.o

$(OBJDIR)/modFACE_parser.o:  $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_help.o $(OBJDIR)/modFACE_misc.o $(OBJDIR)/modFACE_output.o

$(OBJDIR)/modFACE_run.o:  $(OBJDIR)/modFACE_coupling.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_init.o $(OBJDIR)/modFACE_input.o

$(OBJDIR)/modFACE_solver.o:  $(OBJDIR)/modFACE_compute.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_functions.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_output.o $(OBJDIR)/modFACE_precision.o $(OBJDIR)/modFACE_solverf90.o

$(OBJDIR)/modFACE_solverf90.o:  $(OBJDIR)/modFACE_compute.o $(OBJDIR)/modFACE_error.o $(OBJDIR)/modFACE_functions.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_maths.o $(OBJDIR)/modFACE_output.o $(OBJDIR)/modFACE_precision.o

$(OBJDIR)/modFACE_step.o:  $(OBJDIR)/modFACE_compute.o $(OBJDIR)/modFACE_IO.o $(OBJDIR)/modFACE_functions.o $(OBJDIR)/modFACE_header.o $(OBJDIR)/modFACE_output.o $(OBJDIR)/modFACE_solver.o 


