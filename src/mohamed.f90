module PersonModule
    implicit none
    type :: Person
        character(len=100) :: name
        integer :: age
    contains
        procedure :: Say_Hello
    end type Person

    contains
subroutine Say_Hello(this)
    class(Person), intent(in) :: this
    print *, "Hello, my name is ", trim(this%name), " and I am ", this%age, " years old."
end subroutine Say_Hello
end module PersonModule

program main
    use PersonModule
    implicit none
type(Person) :: p 
    p%name = "Mohamed"
    p%age = 30
    call p%Say_Hello()
end program main
