      module modFACE_step
      use modFACE_solver
      use modFACE_header
      use modFACE_functions
      use modFACE_output
      use modFACE_compute

      implicit none
      contains


          subroutine step()

              call shift_array()

              time=time+dt_face
              delta=1.d0/(1.d0+2.d0*pi*nucut*dt_face)

              call compute_source

              call compute_temperature

              call newton_solver

              call check_positivity_max

              call compute_trace_flux

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
