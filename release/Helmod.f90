program Helmod
    use libfemocs, only : femocs
    use iso_c_binding
    implicit none
    
    integer(c_int), parameter :: n_atoms = 1
    real(c_double), dimension(n_atoms) :: x, y, z, phi
    real(c_double), dimension(100000) :: Ex, Ey, Ez, Enorm
    integer(c_int), dimension(1) :: flag = -1
    type(femocs) :: fem
    integer :: counter
    real(c_double) :: t0, t1
    
    do counter = 1, n_atoms
        x(counter) = 48.5 + 1.0 * counter
        y(counter) = 48.5 + 1.0 * counter
        z(counter) = 40.0 + 1.0 * counter
    enddo

    ! Measure the execution time
    call cpu_time(t0)
    call cpu_time(t1)
    
    ! Create the femocs object
    fem = femocs("")

    ! Import the atoms to femocs
    call fem%import_file("")

    ! Run Laplace solver
    call fem%run(0.18d0, "")
    
    ! Export electric field on atoms
    call fem%export_elfield(100000, Ex, Ey, Ez, Enorm)
    
    ! Interpolate electric potential on point with coordinates x,y,z
    call fem%interpolate_phi(n_atoms, x, y, z, phi, flag)
    
    !write(*,*) "Running femocs_speaker..."
    !call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call fem%delete

end program
