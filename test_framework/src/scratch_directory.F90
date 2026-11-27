!> A directory of its own for the files a test case writes and reads.
!|
!| A case that tests file I/O creates one with `create_scratch_directory`, names its files with
!| `scratch_file` (which records them), writes their exact bytes with `write_scratch_file`, and
!| ends with `remove_scratch_directory`, which deletes the recorded files and then the directory.
!| The directory is created in `$TMPDIR` (`/tmp` when unset), named after the process and a
!| counter, so cases never collide and never leave files in the working directory. Removing it
!| fails an assertion when anything else is left in it -- a file the procedure under test wrote
!| that the case did not expect.
!|
!| POSIX `mkdir`, `rmdir` and `getpid` are called through `bind(C)`; Fortran has no directory
!| operations of its own.
!|
!| `check_error_unit_connected` is for a case where a procedure under test fails to open a file:
!| it fails unless Fortran's error unit is still connected afterwards, and connects it again if
!| not, so the rest of the run still reports. `clear_stack` before such a call makes an undefined
!| local unit number read as 0, the error unit, in an unoptimized build, so that closing a unit
!| never opened shows instead of depending on whatever the stack held.
module scratch_directory
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: int32, error_unit
    use asserts, only: assert_true, assert_equal_int
    implicit none
    private
    public :: scratch_directory_type, create_scratch_directory, scratch_file, write_scratch_file
    public :: remove_scratch_directory, clear_stack, check_error_unit_connected

    !> How many files one scratch directory records.
    integer(int32), parameter :: MAX_FILES = 64
    !> The longest path a recorded file may have.
    integer(int32), parameter :: MAX_PATH = 1024

    !> A scratch directory and the files recorded in it.
    type :: scratch_directory_type
        character(len=:), allocatable :: path
            !! The directory, without a trailing slash
        integer(int32) :: n_files = 0
            !! How many of `files` are recorded
        character(len=MAX_PATH), allocatable :: files(:)
            !! The recorded files, as full paths; room for MAX_FILES
    end type scratch_directory_type

    !> How many directories this process has created, to name the next one.
    integer(int32) :: n_created = 0

    interface
        function c_mkdir(path, mode) result(status) bind(C, name="mkdir")
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int) :: status
        end function c_mkdir

        function c_rmdir(path) result(status) bind(C, name="rmdir")
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_rmdir

        function c_getpid() result(pid) bind(C, name="getpid")
            import :: c_int
            integer(c_int) :: pid
        end function c_getpid
    end interface

contains

    !> Creates a new, empty directory for `label` (part of its name, so a leftover shows whose it
    !| was). Fails an assertion if it cannot be created.
    subroutine create_scratch_directory(scratch, label)
        type(scratch_directory_type), intent(out) :: scratch
        character(*), intent(in) :: label
        character(len=MAX_PATH) :: temporary_root, name
        integer(int32) :: root_length, status

        call get_environment_variable("TMPDIR", temporary_root, root_length, status)
        if (status /= 0 .or. root_length == 0) temporary_root = "/tmp"

        allocate (scratch%files(MAX_FILES))
        n_created = n_created + 1
        write (name, '("tox_test_", i0, "_", i0, "_")') c_getpid(), n_created
        scratch%path = trim(temporary_root)//"/"//trim(name)//label
        ! 448 is octal 700: only the owner may enter
        call assert_equal_int(c_mkdir(scratch%path//c_null_char, 448_c_int), 0_int32, &
                              "create_scratch_directory: mkdir "//scratch%path)
    end subroutine create_scratch_directory

    !> The full path of `name` in `scratch`, recorded so that `remove_scratch_directory` deletes it.
    function scratch_file(scratch, name) result(file_path)
        type(scratch_directory_type), intent(inout) :: scratch
        character(*), intent(in) :: name
        character(len=:), allocatable :: file_path

        file_path = scratch%path//"/"//name
        if (any(scratch%files(1:scratch%n_files) == file_path)) return
        call assert_true(scratch%n_files < MAX_FILES .and. len(file_path) <= MAX_PATH, &
                         "scratch_file: room to record "//file_path)
        if (scratch%n_files >= MAX_FILES .or. len(file_path) > MAX_PATH) return
        scratch%n_files = scratch%n_files + 1
        scratch%files(scratch%n_files) = file_path
    end function scratch_file

    !> Writes `content` to `file_path` byte for byte, replacing what was there: no line ending is
    !| added, so a case spells out every `new_line('a')` it wants, and a file without a final line
    !| ending is written as easily as one with. An empty `content` writes an empty file.
    subroutine write_scratch_file(file_path, content)
        character(*), intent(in) :: file_path, content
        integer(int32) :: unit, status

        open (newunit=unit, file=file_path, access='stream', form='unformatted', status='replace', &
              action='write', iostat=status)
        call assert_equal_int(status, 0_int32, "write_scratch_file: open "//file_path)
        if (status /= 0) return
        write (unit, iostat=status) content
        call assert_equal_int(status, 0_int32, "write_scratch_file: write "//file_path)
        close (unit)
    end subroutine write_scratch_file

    !> Deletes every recorded file that exists, then the directory. Fails an assertion if the
    !| directory cannot be removed, which is what happens when something unrecorded is left in it.
    subroutine remove_scratch_directory(scratch)
        type(scratch_directory_type), intent(inout) :: scratch
        integer(int32) :: i_file, unit, status
        logical :: exists

        do i_file = 1, scratch%n_files
            inquire (file=trim(scratch%files(i_file)), exist=exists)
            if (.not. exists) cycle
            open (newunit=unit, file=trim(scratch%files(i_file)), status='old', iostat=status)
            if (status == 0) close (unit, status='delete')
        end do
        scratch%n_files = 0
        call assert_equal_int(c_rmdir(scratch%path//c_null_char), 0_int32, &
                              "remove_scratch_directory: rmdir "//scratch%path//" (is something unexpected left in it?)")
    end subroutine remove_scratch_directory

    !> Overwrites the next 16 KiB of the stack with zeros (see the module's documentation).
    subroutine clear_stack()
        integer(int32), volatile :: filler(4096)

        filler = 0
    end subroutine clear_stack

    !> Fails unless the error unit is connected; connects it to /dev/stderr again if it is not.
    subroutine check_error_unit_connected(msg)
        character(*), intent(in) :: msg
        logical :: connected

        inquire (unit=error_unit, opened=connected)
        call assert_true(connected, msg//": the error unit was closed")
        if (.not. connected) open (unit=error_unit, file="/dev/stderr", action="write")
    end subroutine check_error_unit_connected
end module scratch_directory
