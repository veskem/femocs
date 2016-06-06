program Helmod
    use libfemocs
    use iso_c_binding
    implicit none
    
    integer(c_int), parameter :: n_atoms = 10
    real(c_double), dimension(n_atoms) :: x
    integer(c_int), dimension(n_atoms) :: types
    type(femocs) :: f
    
    x =     (/ 1,2,3,4,5,6,7,8,9,10 /)
    types = (/ 1,1,1,1,1,1,1,1,1,1  /)

    ! Create an object of type foo
    f = femocs("path/to/input/script")

    ! Call bound procedures (member functions)
    write(*,*) "Running femocs.import_atoms"
    call f%import_atoms(0, x, x, x, types)
    
    write(*,*) "Running femocs.run"
    call f%run(10d0)
    
    !write(*,*) "Running femocs_speaker..."
    !call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically here, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call f%delete

end program
