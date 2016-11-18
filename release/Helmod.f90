program Helmod
    use libfemocs, only : femocs
    use iso_c_binding
    implicit none
    
    integer(c_int) :: success
    integer(c_int), parameter :: n_atoms = 1
    real(c_double), dimension(n_atoms) :: x, y, z, phi
    real(c_double), dimension(100000) :: Ex, Ey, Ez, Enorm
    integer(c_int), dimension(1) :: flag = -1
    type(femocs) :: fem
    integer :: counter
    real(c_double) :: t0, t1, coarse_factor
    integer(c_int) :: n_bins
    
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

    ! Import the atoms to femocs
    call fem%import_file(success, "")
    write(*,*) "Result of import_file:", success

    ! Run Laplace solver
    call fem%run(success, 0.18d0, "")
    write(*,*) "Result of run:", success
    
    ! Export electric field on atoms
    call fem%export_elfield(success, 46053, Ex, Ey, Ez, Enorm)
    write(*,*) "Result of export_elfield:", success
    
    ! Interpolate electric potential on point with coordinates x,y,z
    call fem%interpolate_phi(success, n_atoms, x, y, z, phi, flag)
    write(*,*) "Result of interpolate_phi:", success
    
    ! Read command argument from input script
    call fem%parse_int(success, trim("n_bins"), n_bins)
    write(*,*) "Result of parse_int:", success, n_bins
    
    call fem%parse_double(success, trim("coarse_factor"), coarse_factor)
    write(*,*) "Result of parse_double:", success, coarse_factor
    
    !write(*,*) "Running femocs_speaker..."
    !call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call fem%delete

end program
