! Test program to demonstrate the usage of FEMOCS in Fortran
! 2016 Mihkel Veske, University of Helsinki, University of Tartu

program Tester
    use libfemocs, only : femocs
    use iso_c_binding
    implicit none
    
    integer(c_int) :: success, tot_success
    integer(c_int), parameter :: n_atoms = 1
    real(c_double), dimension(n_atoms) :: x, y, z, phi
    real(c_double), dimension(10000) :: Ex, Ey, Ez, Enorm
    integer(c_int), dimension(1) :: flag = -1
    type(femocs) :: fem
    integer :: counter
    real(c_double) :: t0, t1, cmdarg = -1
    
    do counter = 1, n_atoms
        x(counter) = 48.5 + 1.0 * counter
        y(counter) = 48.5 + 1.0 * counter
        z(counter) = 40.0 + 1.0 * counter
    enddo

    ! Measure the execution time
    call cpu_time(t0)
    call cpu_time(t1)
    
    ! Create the femocs object
    fem = femocs("input/md.in")
    tot_success = 0
    
    ! Import the atoms to femocs
    call fem%import_file(success, "")
    tot_success = tot_success + success

    ! Run Laplace solver
    call fem%run(success, 0.2d0, "")
    tot_success = tot_success + success
    
    ! Export electric field on atoms
    call fem%export_elfield(success, 1000, Ex, Ey, Ez, Enorm)
    tot_success = tot_success + success
    
    ! Interpolate electric potential on point with coordinates x,y,z
    call fem%interpolate_phi(success, n_atoms, x, y, z, phi, flag)
    tot_success = tot_success + success
    
    ! Read command argument from input script
    call fem%parse_double(success, trim("Smooth_Factor"), cmdarg)
    
    ! Print quick overview about the success of the run
    write(*,*)
    
    if (tot_success == 0) then; write(*,*) "full run of Femocs      passed"
    else; write(*,*) "full run of Femocs      failed"
    endif
    
    if (cmdarg /= -1.0) then; write(*,*) "reading Smooth_Factor   passed"
    else; write(*,*) "reading Smooth_Factor   failed"
    endif
    
    ! The destructor should be called automatically, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call fem%delete

end program
