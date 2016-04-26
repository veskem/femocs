module libfemocs
    use iso_c_binding

    private
    public :: femocs, femocs_speaker
    
    ! C functions declarations
    interface
        function create_femocs_c(str) bind(C, name="create_femocs")
            use iso_c_binding
            implicit none
            type(c_ptr) :: create_femocs_c
            character(len=1, kind=C_CHAR), intent(in) :: str(*)
        end function

        subroutine delete_femocs_c(femocs) bind(C, name="delete_femocs")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: femocs
        end subroutine

        subroutine femocs_run_c(femocs, E_field, phi) bind(C, name="femocs_run")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            real(c_double), intent(in), value :: E_field
            real(c_double) :: phi(:,:,:)
        end subroutine
        
        subroutine femocs_import_atoms_c(femocs, n_atoms, x, y, z, types) bind(C, name="femocs_import_atoms")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int), value :: n_atoms
            real(c_double) :: x(*)
            real(c_double) :: y(*)
            real(c_double) :: z(*)
            integer(c_int) :: types(*)
        end subroutine

        ! void functions maps to subroutines
        subroutine femocs_speaker_c(str) bind(C, name="femocs_speaker")
            use iso_c_binding
            implicit none
            character(len=1, kind=C_CHAR), intent(in) :: str(*)
        end subroutine
    end interface

    ! We'll use a Fortan type to represent a C++ class here, in an opaque maner
    type femocs
        private
        type(c_ptr) :: ptr ! pointer to the Femocs class
    contains
        ! We can bind some functions to this type, allowing for a cleaner syntax.
        final :: delete_femocs ! Destructor
        procedure :: delete => delete_femocs_polymorph ! Destructor for gfortran
        ! Function members
        procedure :: run => femocs_run
        procedure :: import_atoms => femocs_import_atoms
    end type

    ! This function will act as the constructor for femocs type
    interface femocs
        procedure create_femocs
    end interface

contains ! Implementation of the functions. We just wrap the C function here.
    function create_femocs(str)
        implicit none
        type(femocs) :: create_femocs
        character(len=*), intent(in) :: str
        character(len=1, kind=C_CHAR) :: c_str(len_trim(str) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(str)
        do i = 1, N
            c_str(i) = str(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        create_femocs%ptr = create_femocs_c(c_str)
    end function

    subroutine delete_femocs(this)
        implicit none
        type(femocs) :: this
        call delete_femocs_c(this%ptr)
    end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_femocs_polymorph(this)
        implicit none
        class(femocs) :: this
        call delete_femocs_c(this%ptr)
    end subroutine

    subroutine femocs_run(this, E_field, phi)
        implicit none
        class(femocs), intent(in) :: this
        real(c_double), intent(in) :: E_field
        real(c_double) :: phi(:,:,:)
        call femocs_run_c(this%ptr, E_field, phi)
    end subroutine
    
    subroutine femocs_import_atoms(this, n_atoms, x, y, z, types)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: n_atoms
        real(c_double) :: x(*)
        real(c_double) :: y(*)
        real(c_double) :: z(*)
        integer(c_int) :: types(*)
        call femocs_import_atoms_c(this%ptr, n_atoms, x, y, z, types)
    end subroutine

    subroutine femocs_speaker(str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=1, kind=C_CHAR) :: c_str(len_trim(str) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(str)
        do i = 1, N
            c_str(i) = str(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call femocs_speaker_c(c_str)
    end subroutine
end module
