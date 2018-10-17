      module modFACE_step
      use modFACE_solver
      use modFACE_header
      use modFACE_functions
      use modFACE_output
      use modFACE_compute
      use modFACE_IO
      implicit none
      contains


          subroutine step()
              integer::k

              logical :: quick_convergence

              if (verbose_debug) call print_milestone('step')
              ! update density and temp ndt->ndt-1
              call shift_array
                do k=1,nspc
            if (inflx_in_pulse(k).ne."N") then
                if ((time+dt_face.ge.inflx_in_pulse_starttime(k)).and.(time.lt.inflx_in_pulse_starttime(k))) then
                dt_face=inflx_in_pulse_starttime(k)-time
                endif
                endif
                enddo

              ! calculate new source and temperature both because of external influx with time dependecy)
              call compute_source
              call compute_temperature


              ! ***** now  solve dndt=f(n(t+dt)) equation *****

              quick_convergence=.true.
                !  quick_convergence  F : no quick convergence and no update of density and temp in newton_solver routine
                !                     T : solver step converged and update of density and temp in newton_solver routine



              ! loop to adjust dt_face ot obtain fast inversion of jacobian
              !  while solver convergence is not satisfying (solver_stauts.ne.1) and dt_face is not equal to min_dt_face
  !            do while(.not.quick_convergence)
   ! if reduction factor =1 then no time step reduction
!               if (reduction_factor_dt.eq.1d0) then
                   !   quick_convergence=.true.
!               endif

                  ! get delta for low-freq filter
                  !delta=1.d0/(1.d0+2.d0*pi*nucut*dt_face)


                      ! if dt_face is smaller than the minimun timestep requested then
                      ! use minimal timestep and force update of density and temp after step (quick_convergence=T)
!                      if (dt_face.le.min_dt_face) then
!                          dt_face=min_dt_face
!                          quick_convergence=.true.
!                      endif

                  call newton_solver(quick_convergence)

                  ! if solver struggles to get fast converging solution then reduction of the timestep for the next iteration of the while loop

!                  if (.not.quick_convergence) then
!                      dt_face=dt_face*reduction_factor_dt
!                      call print_reduction_timestep_info
!                 endif

!              enddo


              ! we count the number of successful steps (without time step reduction)
              solver_step_count=solver_step_count+1

              ! update time when solver step is completed with the current dt_face
              time=time+dt_face
              time_savetime=time_savetime+dt_face
              time_savevol=time_savevol+dt_face
              call check_positivity_max
              call compute_trace_flux
              call compute_onthefly_inventory
              call print_timestep_info                   ! print info on current time step
              call compute_dt_update
                            ! if dt_face is less than nominal time step dt0_face,
                           ! then we increase dt by 1/reduction_factor_dt if there was enough successsful timestep
                           ! for the next iteration
!              if (dt_face.lt.dt0_face.and.solver_step_count.gt.Nstep_increase_dt.and.reduction_factor_dt.ne.1d0) then
!                  solver_step_count=0 ! we restart the count to let the solver run several time
!                  dt_face=dt_face/reduction_factor_dt ! we increase the time step
!                  if (dt_face.gt.dt0_face) dt_face=dt0_face ! we
!              endif
          end subroutine step


    subroutine shift_array
    ! shifting time array down in time (current time: ndt, past time: ndt-1,ndt-2,...)
        integer i,j,k,l
        ! volume
        do i=1,ndt-1
            do j=0,ngrd
                temp(i,j)=temp(i+1,j)
                rate_t (i,j)=rate_t (i+1,j)
                qflx(i,j)=qflx(i+1,j)
                do k=1,nspc
                    dens(i,j,k)=dens(i+1,j,k)
                    flx (i,j,k)=flx (i+1,j,k)
                    ero_flx (i,j,k)=ero_flx (i+1,j,k)
                    dif_flx (i,j,k)=dif_flx (i+1,j,k)
                    src (i,j,k)=src (i+1,j,k)
                    srs (i,j,k)=srs (i+1,j,k)
                    cdif(i,j,k)=cdif(i+1,j,k)
                    rct (i,j,k)=rct (i+1,j,k)
                    rate_d (i,j,k)=rate_d (i+1,j,k)
                    do l=1,nspc
                        srb(i,j,k,l)=srb(i+1,j,k,l)
                    enddo
                enddo
            enddo
        enddo
        ! surface
        do i=1,ndt-1
            do k=1,nspc
                dsrfl(i,k)=dsrfl(i+1,k)
                Gsrf_l (i,k)=Gsrf_l (i+1,k)
                Gabs_l  (i,k)=Gabs_l  (i+1,k)
                Gdes_l  (i,k)=Gdes_l  (i+1,k)
                Gb_l  (i,k)=Gb_l  (i+1,k)
                Gads_l  (i,k)=Gads_l  (i+1,k)
                dsrfr(i,k)=dsrfr(i+1,k)
                Gsrf_r (i,k)=Gsrf_r (i+1,k)
                Gabs_r  (i,k)=Gabs_r  (i+1,k)
                Gdes_r  (i,k)=Gdes_r  (i+1,k)
                Gb_r  (i,k)=Gb_r  (i+1,k)
                Gads_r  (i,k)=Gads_r  (i+1,k)
                jout (i,k)=jout (i+1,k)
            enddo
        enddo

    end subroutine shift_array





      end module
