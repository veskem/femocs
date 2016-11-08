program Helmod
    use libfemocs, only : femocs
    use iso_c_binding
    implicit none
    
    integer(c_int), parameter :: n_atoms = 10
    real(c_double), dimension(n_atoms) :: x
    integer(c_int), dimension(n_atoms) :: types
    type(femocs) :: fem
    
    x =     (/ 1,2,3,4,5,6,7,8,9,10 /)
    types = (/ 1,1,1,1,1,1,1,1,1,1  /)

    ! Create an object of type foo
    fem = femocs("")

    ! Call bound procedures (member functions)
    write(*,*) "Running femocs.import_atoms"
    ! call fem%import_atoms(0, x, x, x, types)
    call fem%import_file("")
     
    write(*,*) "Running femocs.run"
    call fem%run(0.1d0, "")
    
    !write(*,*) "Running femocs_speaker..."
    !call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically here, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call fem%delete

end program
