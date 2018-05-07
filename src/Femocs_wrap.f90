module libfemocs
    use iso_c_binding

    private
    public :: femocs
    
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

        subroutine femocs_run_c(femocs, retval, timestep) bind(C, name="femocs_run")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval
            integer(c_int), intent(in), value :: timestep
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
        
        subroutine femocs_interpolate_surface_elfield_c(femocs,retval,n_points,x,y,z,Ex,Ey,Ez,Enorm,flag) &
                                                 bind(C, name="femocs_interpolate_surface_elfield")
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

        subroutine femocs_interpolate_phi_c(femocs,retval,n_points,x,y,z,phi,flag) &
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

        subroutine femocs_export_results_c(femocs,retval,n_points,data_type,data) &
                                                 bind(C, name="femocs_export_results")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval
            integer(c_int), value :: n_points
            character(len=1, kind=C_CHAR), intent(in) :: data_type(*)
            real(c_double) :: data(*)
        end subroutine

        subroutine femocs_interpolate_results_c(femocs,retval,n_points,data_type,near_surface,x,y,z,data,flag) &
                                                    bind(C, name="femocs_interpolate_results")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: femocs
            integer(c_int) :: retval
            integer(c_int), value :: n_points
            character(len=1, kind=C_CHAR), intent(in) :: data_type(*)
            integer(c_int), value :: near_surface
            real(c_double) :: x(*)
            real(c_double) :: y(*)
            real(c_double) :: z(*)
            real(c_double) :: data(*)
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
        ! Destructor for gfortran
        procedure :: delete => delete_femocs
        ! Function members
        procedure :: run => femocs_run
        procedure :: import_atoms => femocs_import_atoms
        procedure :: import_parcas => femocs_import_parcas
        procedure :: import_file => femocs_import_file
        procedure :: export_atom_types => femocs_export_atom_types
        procedure :: interpolate_elfield => femocs_interpolate_elfield
        procedure :: interpolate_surface_elfield => femocs_interpolate_surface_elfield
        procedure :: interpolate_phi => femocs_interpolate_phi
        procedure :: export_results => femocs_export_results
        procedure :: interpolate_results => femocs_interpolate_results
        procedure :: parse_int => femocs_parse_int
        procedure :: parse_double => femocs_parse_double
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

    subroutine delete_femocs(this)
        implicit none
        class(femocs) :: this
        call delete_femocs_c(this%ptr)
    end subroutine

    subroutine femocs_run(this, retval, timestep)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval
        integer(c_int), intent(in) :: timestep
        call femocs_run_c(this%ptr, retval, timestep)
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

    subroutine femocs_interpolate_surface_elfield(this,retval,n_points,x,y,z,Ex,Ey,Ez,Enorm,flag)
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
        call femocs_interpolate_surface_elfield_c(this%ptr,retval,n_points,x,y,z,Ex,Ey,Ez,Enorm,flag)
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

    subroutine femocs_export_results(this, retval, n_points, data_type, data)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval
        integer(c_int) :: n_points
        character(len=*), intent(in) :: data_type
        real(c_double) :: data(*)
        character(len=1, kind=C_CHAR) :: c_str(len_trim(data_type) + 1)
        integer :: N, i

        ! Convert Fortran string to C string
        N = len_trim(data_type)
        do i = 1, N
            c_str(i) = data_type(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call femocs_export_results_c(this%ptr, retval, n_points, c_str,data)
    end subroutine

    subroutine femocs_interpolate_results(this,retval,n_points,data_type,near_surface,x,y,z,data,flag)
        implicit none
        class(femocs), intent(in) :: this
        integer(c_int) :: retval
        integer(c_int) :: n_points
        character(len=*), intent(in) :: data_type
        integer(c_int) :: near_surface
        real(c_double) :: x(*)
        real(c_double) :: y(*)
        real(c_double) :: z(*)
        real(c_double) :: data(*)
        integer(c_int) :: flag(*)
        character(len=1, kind=C_CHAR) :: c_str(len_trim(data_type) + 1)
        integer :: N, i

        ! Convert Fortran string to C string
        N = len_trim(data_type)
        do i = 1, N
            c_str(i) = data_type(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call femocs_interpolate_results_c(this%ptr,retval,n_points,c_str,near_surface,x,y,z,data,flag)
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
