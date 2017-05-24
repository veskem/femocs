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

        subroutine femocs_run_c(femocs, retval, E_field, message) bind(C, name="femocs_run")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval
            real(c_double), intent(in), value :: E_field
            character(len=1, kind=C_CHAR), intent(in) :: message(*)
        end subroutine
       
        subroutine femocs_import_atoms_c(femocs, retval, n_atoms, x, y, z, types) bind(C, name="femocs_import_atoms")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval
            integer(c_int), value :: n_atoms
            real(c_double) :: x(*)
            real(c_double) :: y(*)
            real(c_double) :: z(*)
            integer(c_int) :: types(*)
        end subroutine

        subroutine femocs_import_parcas_c(femocs, retval, n_atoms, coordinates, box, nborlist) bind(C, name="femocs_import_parcas")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_atoms
            real(c_double) :: coordinates(*)
            real(c_double) :: box(*)
            integer(c_int) :: nborlist(*)
        end subroutine
        
        subroutine femocs_import_file_c(femocs, retval, file_name) bind(C, name="femocs_import_file")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            character(len=1, kind=C_CHAR), intent(in) :: file_name(*)
        end subroutine
        
        subroutine femocs_export_atom_types_c(femocs, retval, n_atoms, types) bind(C, name="femocs_export_atom_types")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_atoms
            integer(c_int) :: types(*)
        end subroutine
        
        subroutine femocs_export_elfield_c(femocs, retval, n_atoms, Ex, Ey, Ez, Enorm) bind(C, name="femocs_export_elfield")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_atoms
            real(c_double) :: Ex(*)
            real(c_double) :: Ey(*)
            real(c_double) :: Ez(*)
            real(c_double) :: Enorm(*)
        end subroutine
        
        subroutine femocs_export_temperature_c(femocs, retval, n_atoms, T) bind(C, name="femocs_export_temperature")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_atoms
            real(c_double) :: T(*)
        end subroutine
        
        subroutine femocs_export_charge_and_force_c(femocs, retval, n_atoms, xq) bind(C, name="femocs_export_charge_and_force")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_atoms
            real(c_double) :: xq(*)
        end subroutine
        
        subroutine femocs_interpolate_elfield_c(femocs, retval, n_points, x, y, z, Ex, Ey, Ez, Enorm, flag) &
                                                 bind(C, name="femocs_interpolate_elfield")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_points
            real(c_double) :: x(*)
            real(c_double) :: y(*)
            real(c_double) :: z(*)
            real(c_double) :: Ex(*)
            real(c_double) :: Ey(*)
            real(c_double) :: Ez(*)
            real(c_double) :: Enorm(*)
            integer(c_int) :: flag(*)
        end subroutine
        
        subroutine femocs_interpolate_phi_c(femocs, retval, n_points, x, y, z, phi, flag) &
                                                 bind(C, name="femocs_interpolate_phi")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            integer(c_int), value :: n_points
            real(c_double) :: x(*)
            real(c_double) :: y(*)
            real(c_double) :: z(*)
            real(c_double) :: phi(*)
            integer(c_int) :: flag(*)
        end subroutine
        
        subroutine femocs_parse_int_c(femocs, retval, command, arg) bind(C, name="femocs_parse_int")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            character(len=1, kind=C_CHAR), intent(in) :: command(*)
            integer(c_int) :: arg
        end subroutine
        
        subroutine femocs_parse_double_c(femocs, retval, command, arg) bind(C, name="femocs_parse_double")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            character(len=1, kind=C_CHAR), intent(in) :: command(*)
            real(c_double) :: arg
        end subroutine
        
        subroutine femocs_parse_boolean_c(femocs, retval, command, arg) bind(C, name="femocs_parse_boolean")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            character(len=1, kind=C_CHAR), intent(in) :: command(*)
            logical(c_bool) :: arg
        end subroutine
        
        subroutine femocs_parse_string_c(femocs, retval, command, arg) bind(C, name="femocs_parse_string")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval            
            character(len=1, kind=C_CHAR), intent(in) :: command(*)
            character(len=1, kind=C_CHAR) :: arg(*)
        end subroutine
        
    end interface

    ! We'll use a Fortan type to represent a C++ class here in an opaque manner
    type femocs
        private
        type(c_ptr) :: ptr ! pointer to the Femocs class
    contains
        ! Compiler gives the following warning when this line is uncommented: 
        ! Only array FINAL procedures declared for derived type ‘femocs’ defined at (1), suggest also scalar one [-Wsurprising]
        ! For some reason the original sample of the wrapper had that function and I'm not exactly sure what it does, therefore I didn't dare to erase it.
        ! However, it seems wrapper can perfectly do without it
!         final :: delete_femocs

        ! Destructor for gfortran
        procedure :: delete => delete_femocs_polymorph
        ! Function members
        procedure :: run => femocs_run
        procedure :: import_atoms => femocs_import_atoms
        procedure :: import_parcas => femocs_import_parcas
        procedure :: import_file => femocs_import_file
        procedure :: export_atom_types => femocs_export_atom_types
        procedure :: export_elfield => femocs_export_elfield
        procedure :: export_temperature => femocs_export_temperature
        procedure :: export_charge_and_force => femocs_export_charge_and_force
        procedure :: interpolate_elfield => femocs_interpolate_elfield
        procedure :: interpolate_phi => femocs_interpolate_phi
        procedure :: parse_int => femocs_parse_int
        procedure :: parse_double => femocs_parse_double
        procedure :: parse_boolean => femocs_parse_boolean
        procedure :: parse_string => femocs_parse_string
    end type       
        
    ! This function will act as the constructor for femocs type
    interface femocs
        procedure create_femocs
    end interface

    ! Implementation of the functions. We just wrap the C function here.
    contains 
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

    ! See the comment about compiler warnings above
!     subroutine delete_femocs(this)
!         implicit none
!         type(femocs) :: this
!         call delete_femocs_c(this%ptr)
!     end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_femocs_polymorph(this)
        implicit none
        class(femocs) :: this
        call delete_femocs_c(this%ptr)
    end subroutine

    subroutine femocs_run(this, retval, E_field, message)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval
        real(c_double), intent(in) :: E_field
        character(len=*), intent(in) :: message
        character(len=1, kind=C_CHAR) :: c_str(len_trim(message) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(message)
        do i = 1, N
            c_str(i) = message(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call femocs_run_c(this%ptr, retval, E_field, c_str)
    end subroutine
    
    subroutine femocs_import_atoms(this, retval, n_atoms, x, y, z, types)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval
        integer(c_int) :: n_atoms
        real(c_double) :: x(*)
        real(c_double) :: y(*)
        real(c_double) :: z(*)
        integer(c_int) :: types(*)
        call femocs_import_atoms_c(this%ptr, retval, n_atoms, x, y, z, types)
    end subroutine
    
    subroutine femocs_import_parcas(this, retval, n_atoms, coordinates, box, nborlist)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        integer(c_int) :: n_atoms
        real(c_double) :: coordinates(*)
        real(c_double) :: box(*)
        integer(c_int) :: nborlist(*)
        call femocs_import_parcas_c(this%ptr, retval, n_atoms, coordinates, box, nborlist)
    end subroutine    

    subroutine femocs_import_file(this, retval, file_name)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        character(len=*), intent(in) :: file_name
        character(len=1, kind=C_CHAR) :: c_str(len_trim(file_name) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(file_name)
        do i = 1, N
            c_str(i) = file_name(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call femocs_import_file_c(this%ptr, retval, c_str)
    end subroutine
    
    subroutine femocs_export_atom_types(this, retval, n_atoms, types)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        integer(c_int) :: n_atoms
        integer(c_int) :: types(*)
        call femocs_export_atom_types_c(this%ptr, retval, n_atoms, types)
    end subroutine
    
    subroutine femocs_export_elfield(this, retval, n_atoms, Ex, Ey, Ez, Enorm)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        integer(c_int) :: n_atoms
        real(c_double) :: Ex(*)
        real(c_double) :: Ey(*)
        real(c_double) :: Ez(*)
        real(c_double) :: Enorm(*)
        call femocs_export_elfield_c(this%ptr, retval, n_atoms, Ex, Ey, Ez, Enorm)
    end subroutine

    subroutine femocs_export_temperature(this, retval, n_atoms, T)
        implicit none
        class(femocs), intent(in) :: this
        type(c_ptr), intent(in), value :: femocs
        integer(c_int) :: retval            
        integer(c_int), value :: n_atoms
        real(c_double) :: T(*)
        call femocs_export_temperature_c(this%ptr, retval, n_atoms, T)
    end subroutine
    
    subroutine femocs_export_charge_and_force(this, retval, n_atoms, xq)
        implicit none
        class(femocs), intent(in) :: this
        type(c_ptr), intent(in), value :: femocs
        integer(c_int) :: retval            
        integer(c_int), value :: n_atoms
        real(c_double) :: xq(*)
        call femocs_export_charge_and_force_c(this%ptr, retval, n_atoms, xq)
    end subroutine
        
    subroutine femocs_interpolate_elfield(this, retval, n_points, x, y, z, Ex, Ey, Ez, Enorm, flag)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        integer(c_int) :: n_points
        real(c_double) :: x(*)
        real(c_double) :: y(*)
        real(c_double) :: z(*)
        real(c_double) :: Ex(*)
        real(c_double) :: Ey(*)
        real(c_double) :: Ez(*)
        real(c_double) :: Enorm(*)
        integer(c_int) :: flag(*)
        call femocs_interpolate_elfield_c(this%ptr, retval, n_points, x, y, z, Ex, Ey, Ez, Enorm, flag)
    end subroutine
    
    subroutine femocs_interpolate_phi(this, retval, n_points, x, y, z, phi, flag)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        integer(c_int) :: n_points
        real(c_double) :: x(*)
        real(c_double) :: y(*)
        real(c_double) :: z(*)
        real(c_double) :: phi(*)
        integer(c_int) :: flag(*)
        call femocs_interpolate_phi_c(this%ptr, retval, n_points, x, y, z, phi, flag)
    end subroutine

    subroutine femocs_parse_int(this, retval, command, arg)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        character(len=*), intent(in) :: command
        integer(c_int) :: arg
        character(len=1, kind=C_CHAR) :: c_str(len_trim(command) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(command)
        do i = 1, N
            c_str(i) = command(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR       
        
        call femocs_parse_int_c(this%ptr, retval, c_str, arg)
    end subroutine
        
    subroutine femocs_parse_double(this, retval, command, arg)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        character(len=*), intent(in) :: command
        real(c_double) :: arg
        character(len=1, kind=C_CHAR) :: c_str(len_trim(command) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(command)
        do i = 1, N
            c_str(i) = command(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR       
        
        call femocs_parse_double_c(this%ptr, retval, c_str, arg)
    end subroutine
    
    subroutine femocs_parse_boolean(this, retval, command, arg)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        character(len=*), intent(in) :: command
        logical(c_bool) :: arg
        character(len=1, kind=C_CHAR) :: c_str(len_trim(command) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(command)
        do i = 1, N
            c_str(i) = command(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR       
        
        call femocs_parse_boolean_c(this%ptr, retval, c_str, arg)
    end subroutine
    
    subroutine femocs_parse_string(this, retval, command, arg)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval        
        character(len=*), intent(in) :: command
        character(len=*) :: arg
        character(len=1, kind=C_CHAR) :: c_str(len_trim(command) + 1)
        character(len=1, kind=C_CHAR) :: c_arg(len_trim(arg) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(command)
        do i = 1, N
            c_str(i) = command(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR       
        
        call femocs_parse_string_c(this%ptr, retval, c_str, c_arg)
        
        ! Converting C argument to Fortran argument
        N = len_trim(arg)
        do i = 1, N
            arg(i:i) = c_arg(i)
        end do   
    end subroutine    

end module
