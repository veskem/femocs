
! Test program to call Femocs from Fortran

program main
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine femocs_speaker() bind (c)
      use iso_c_binding
    end subroutine femocs_speaker
  end interface

  write ( *, '(a)' ) '=== Helmod calling Femocs ==='
  write ( *, '(a)' ) ' '
    
  call femocs_speaker()

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '=== Helmod ended ==='

  stop
end
