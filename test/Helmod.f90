program Helmod2
    use libfemocs
    use iso_c_binding
    implicit none
    
    integer(c_int), parameter :: n_atoms = 10
    integer :: i, j, k
    real(c_double), dimension(n_atoms) :: x
    real(c_double), dimension(n_atoms,n_atoms,n_atoms) :: phi
    type(femocs) :: f
    
    x = (/ 1,2,3,4,5,6,7,8,9,10 /)
    
    do i = 1, n_atoms
      do j = 1, n_atoms
        do k = 1, n_atoms
          phi(i,j,k) = i*j*k
        end do
      end do
    end do

    ! Create an object of type foo
    f = femocs("path/to/input/script")

    ! Call bound procedures (member functions)
    write(*,*) "Running femocs.import_atoms"
    call f%import_atoms(0, x, x, x)
    
    write(*,*) "Running femocs.run"
    call f%run(10d0, phi)
    
!    do i = 1, n_atoms
!      do j = 1, n_atoms
!        do k = 1, n_atoms
!          write(*,*) 'phi',i,j,k, phi(i,j,k)
!        end do
!      end do
!    end do
    
    !write(*,*) "Running femocs_speaker..."
    !call femocs_speaker("From Fortran!")

    ! The destructor should be called automatically here, but this is not yet
    ! implemented in gfortran. So let's do it manually.
    call f%delete

end program
