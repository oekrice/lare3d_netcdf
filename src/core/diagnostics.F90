!  Copyright 2020 University of Warwick

!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at

!      http://www.apache.org/licenses/LICENSE-2.0

!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.

!******************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare3d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE version_data
  USE netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_diags, output_snap, energy_correction

  REAL(dbl) :: visc_heating
  LOGICAL, SAVE :: visc_heating_updated = .FALSE.

CONTAINS

  !****************************************************************************
  ! My diagnostic functions (same as for mf2d)
  !****************************************************************************

  SUBROUTINE output_diags(diag_num)
      !Calculates some diagnostics and saves to netcdf file as for the triangle code, which was fairly neat (if i say so myself...). Should make for easy pythonning.
      IMPLICIT NONE
      INTEGER:: diag_num
      REAL(num):: diag
      character(len=100):: filename
      integer:: id_1, id_2, id_3, id_4, id_5, id_6, id_7, id_8, ncid, nd_id

      !Allocate diagnostic arrays
      if (diag_num == 0) then
          allocate(diag_time(0:ndiags))
      end if

      !TIME
      diag = time
      diag_time(diag_num) = time

      !ADMIN
      if (run_id < 10) then
      write (filename, "(A19,I1,A3)") "./diagnostics/run00", int(run_id), ".nc"
      else if (run_id < 100) then
          write (filename, "(A18,I2,A3)") "./diagnostics/run0", int(run_id), ".nc"
      else
          write (filename, "(A17,I3,A3)") "./diagnostics/run", int(run_id), ".nc"
      end if

      !Write to diagnostics file, using netcdf
      call try(nf90_create(trim(filename), nf90_clobber, ncid))
      call try(nf90_def_dim(ncid, 'ndiags', ndiags+1, nd_id))  !Make up fake dimensions here

      call try(nf90_def_var(ncid, 'time', nf90_double, (/nd_id/), id_1))
      call try(nf90_enddef(ncid))

      call try(nf90_put_var(ncid, id_1, diag_time))


      call try(nf90_close(ncid))

      diag_num = diag_num + 1

  END SUBROUTINE output_diags


SUBROUTINE output_snap(snap_num)
    !Exports the magnetic field at this plot_num to an appropriate netcdf file
    IMPLICIT NONE

    CHARACTER(LEN =64):: output_filename
    INTEGER:: snap_num, proc_write
    INTEGER:: ncid, vid
    INTEGER:: xs_id, ys_id, zs_id
    INTEGER:: xc_id, yc_id, zc_id
    INTEGER:: bx_id, by_id, bz_id
    INTEGER:: en_id, rho_id

    if (snap_num < 10) then
        write (output_filename, "(A7,A3,I1,A3)") "./Data/", "000", snap_num, ".nc"
    else if (snap_num < 100) then
        write (output_filename, "(A7,A2,I2,A3)") "./Data/", "00", snap_num, ".nc"
    else if (snap_num < 1000) then
        write (output_filename, "(A7,A1,I3,A3)")  "./Data/", "0", snap_num, ".nc"
    else if (snap_num < 10000) then
        write (output_filename, "(A7,I4,A3)")  "./Data/", snap_num, ".nc"
    end if

    if (rank == 0) then
    call try(nf90_create(trim(output_filename), nf90_clobber, ncid))

    call try(nf90_def_dim(ncid, 'xs', nx_global+1, xs_id))
    call try(nf90_def_dim(ncid, 'ys', ny_global+1, ys_id))
    call try(nf90_def_dim(ncid, 'zs', nz_global+1, zs_id))

    call try(nf90_def_dim(ncid, 'xc', nx_global, xc_id))
    call try(nf90_def_dim(ncid, 'yc', ny_global, yc_id))
    call try(nf90_def_dim(ncid, 'zc', nz_global, zc_id))

    call try(nf90_def_var(ncid, 'bx', nf90_double, (/xs_id ,yc_id, zc_id/), bx_id))
    call try(nf90_def_var(ncid, 'by', nf90_double, (/xc_id ,ys_id, zc_id/), by_id))
    call try(nf90_def_var(ncid, 'bz', nf90_double, (/xc_id ,yc_id, zs_id/), bz_id))

    call try(nf90_def_var(ncid, 'en', nf90_double, (/xc_id ,yc_id, zc_id/), en_id))
    call try(nf90_def_var(ncid, 'rho', nf90_double, (/xc_id ,yc_id, zc_id/), rho_id))

    call try(nf90_enddef(ncid))
    call try(nf90_close(ncid))

    end if
    call MPI_BARRIER(comm,errcode)

    !Each process writes data in turn

    do proc_write = 0 ,nproc-1
        call MPI_BARRIER(comm,errcode)

        if (rank == proc_write) then
            call try(nf90_open(trim(output_filename), nf90_write, ncid))

            call try(nf90_inq_varid(ncid, 'bx', vid))
            call try(nf90_put_var(ncid, vid, bx(0:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'by', vid))
            call try(nf90_put_var(ncid, vid, by(1:nx,0:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'bz', vid))
            call try(nf90_put_var(ncid, vid, bz(1:nx,1:ny,0:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'en', vid))
            call try(nf90_put_var(ncid, vid, energy(1:nx,1:ny,1:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'rho', vid))
            call try(nf90_put_var(ncid, vid, rho(1:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))

            call try(nf90_close(ncid))

        end if
        call MPI_BARRIER(comm,errcode)

    end do

    call mpi_barrier(comm, errcode)
    if (rank == 0) print*, 'Saved snapshot number', snap_num, ' at time', time, 'to file ', output_filename

    snap_num = snap_num + 1
    return


END SUBROUTINE output_snap


SUBROUTINE try(status)
  ! Catch error in reading netcdf field.
  INTEGER, INTENT(IN):: status

  if (status /= NF90_noerr) THEN
      PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
  end if

END SUBROUTINE try



  !****************************************************************************
  ! Test whether any of the conditions for doing output on the current
  ! iteration are met
  !****************************************************************************

  SUBROUTINE io_test(step, print_arrays, last_call)

    INTEGER, INTENT(IN) :: step
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time_prev + dt_snapshots
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      time_prev = t1
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. step == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

    IF (restart) THEN
      print_arrays = .FALSE.
      file_number = file_number + 1
      restart = .FALSE.
    END IF

  END SUBROUTINE io_test




  SUBROUTINE energy_account(energy_b, energy_ke, energy_int, do_sum)

    REAL(dbl), INTENT(OUT) :: energy_b, energy_ke, energy_int
    LOGICAL, INTENT(IN) :: do_sum
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(dbl) :: energy_local(3), energy_sum(3)
    REAL(dbl) :: cv_v, rho_v, w1, w2, w3

    energy_b_local   = 0.0_dbl
    energy_ke_local  = 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          w1 = (bx(ix,iy,iz)**2 + bx(ixm,iy ,iz )**2) * 0.5_num
          w2 = (by(ix,iy,iz)**2 + by(ix ,iym,iz )**2) * 0.5_num
          w3 = (bz(ix,iy,iz)**2 + bz(ix ,iy ,izm)**2) * 0.5_num
          w1 = (w1 + w2 + w3) * 0.5_dbl
          energy_b_local = energy_b_local + w1 * cv(ix,iy,iz)

          energy_int_local = energy_int_local &
              + energy(ix,iy,iz) * rho(ix,iy,iz) * cv(ix,iy,iz)
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          ! WARNING the KE is summed on the vertices
          rho_v = rho(ix ,iy ,iz ) * cv(ix ,iy ,iz ) &
                + rho(ixp,iy ,iz ) * cv(ixp,iy ,iz ) &
                + rho(ix ,iyp,iz ) * cv(ix ,iyp,iz ) &
                + rho(ixp,iyp,iz ) * cv(ixp,iyp,iz ) &
                + rho(ix ,iy ,izp) * cv(ix ,iy ,izp) &
                + rho(ixp,iy ,izp) * cv(ixp,iy ,izp) &
                + rho(ix ,iyp,izp) * cv(ix ,iyp,izp) &
                + rho(ixp,iyp,izp) * cv(ixp,iyp,izp)

          cv_v = cv(ix,iy ,iz ) + cv(ixp,iy ,iz ) &
               + cv(ix,iyp,iz ) + cv(ixp,iyp,iz ) &
               + cv(ix,iy ,izp) + cv(ixp,iy ,izp) &
               + cv(ix,iyp,izp) + cv(ixp,iyp,izp)

          rho_v = rho_v / cv_v
          cv_v = cv_v * 0.125_dbl
          w1 = rho_v * cv_v &
              * (vx(ix,iy,iz)**2 + vy(ix,iy,iz)**2 + vz(ix,iy,iz)**2)

          IF (ix == 0 .OR. ix == nx) THEN
            w1 = w1 * 0.5_dbl
          END IF

          IF (iy == 0 .OR. iy == ny) THEN
            w1 = w1 * 0.5_dbl
          END IF

          IF (iz == 0 .OR. iz == nz) THEN
            w1 = w1 * 0.5_dbl
          END IF

          energy_ke_local = energy_ke_local + w1 * 0.5_dbl
        END DO
      END DO
    END DO

    IF (do_sum) THEN
      energy_local(1) = energy_b_local
      energy_local(2) = energy_ke_local
      energy_local(3) = energy_int_local

      CALL MPI_ALLREDUCE(energy_local, energy_sum, 3, MPI_DOUBLE_PRECISION, &
          MPI_SUM, comm, errcode)

      energy_b   = energy_sum(1)
      energy_ke  = energy_sum(2)
      energy_int = energy_sum(3)
    ELSE
      energy_b   = energy_b_local
      energy_ke  = energy_ke_local
      energy_int = energy_int_local
    END IF

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke(:,:,:) = delta_ke(:,:,:) / (rho(:,:,:) * cv(:,:,:))

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix,iy,iz) = energy(ix,iy,iz) + delta_ke(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL energy_bcs

  END SUBROUTINE energy_correction

  SUBROUTINE output_log

    ! Writes basic data to 'lare3d.dat'

    IF (restart) THEN
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) '#####################################################'
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) 'Restarting from step ', step, ' and time ', time
      WRITE(stat_unit,*)
    END IF

    WRITE(stat_unit,*) ascii_header
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'nprocx, nprocy, nprocz = ', nprocx, nprocy, nprocz
    WRITE(stat_unit,*) 'nx, ny, nz = ', nx_global, ny_global, nz_global
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'length_x = ', length_x
    WRITE(stat_unit,*) 'length_y = ', length_y
    WRITE(stat_unit,*) 'length_z = ', length_z
    WRITE(stat_unit,*)
#ifdef QMONO
    WRITE(stat_unit,*) 'q_mono viscosity (-DQMONO)'
#else
    WRITE(stat_unit,*) 'tensor shock viscosity'
#endif
#ifdef FOURTHORDER
    WRITE(stat_unit,*) '4th-order resistive update (-DFOURTHORDER)'
#endif
#ifdef SINGLE
    WRITE(stat_unit,*) 'single precision (-DSINGLE)'
#endif
    WRITE(stat_unit,*) 'linear viscosity coeff = ', visc1
    WRITE(stat_unit,*) 'quadratic viscosity coeff = ', visc2
    WRITE(stat_unit,*) 'j_max = ', j_max
    WRITE(stat_unit,*) 'eta0 = ', eta0
    WRITE(stat_unit,*) 'eta_background = ', eta_background
    WRITE(stat_unit,*) 'kappa = ', kappa_0
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'mass_fraction = ', mf
    WRITE(stat_unit,*) 'normalising B = ', B_norm
    WRITE(stat_unit,*) 'normalising L = ', L_norm
    WRITE(stat_unit,*) 'normalising density = ', rho_norm
    WRITE(stat_unit,*) 'normalising speed = ', B_norm / SQRT(mu0_SI * rho_norm)
    WRITE(stat_unit,*) 'normalising time = ', L_norm / (B_norm / SQRT(mu0_SI * rho_norm))
    WRITE(stat_unit,*) 'normalising temperature = ', temp_norm
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 't_start, t_end = ', time, t_end
    WRITE(stat_unit,*) 'nsteps =', nsteps
    WRITE(stat_unit,*)

  END SUBROUTINE output_log

END MODULE diagnostics
