program Helmod2
    use libfemocs
    implicit none
    ! type(femocs) :: f

    ! Create an object of type foo
    ! f = femocs("path/to/input/script")

    ! Call bound procedures (member functions)
    ! write(*,*) "Running femocs.run(10.0)"
    ! call f%run(10d0)
    
    write(*,*) "Running femocs_speaker..."
    call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically here, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    ! call f%delete

end program
