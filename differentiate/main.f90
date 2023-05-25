#define _pr 8

module config
    implicit none
    public

    !Differentiation step
    real(_pr), parameter :: del = .00001
    !Segment boundaries
    real(_pr), parameter :: min = 0, max = 10
    !Plot smoothness
    integer, parameter :: dots = 100
end module

program draw
    use config
    implicit none

    integer :: i
    real(_pr) :: x
    procedure(func), pointer :: fptr
    fptr => func

    open(unit = 1, file = "output.dat")
    do i = 1, dots
        x = min + (i-1)*((max-min)/(dots-1))
        write(1, *) x, func(x), diff(x, fptr), ddiff(x, fptr)
    end do
    close(1)

contains
    function func(arg) result(res)
        real(_pr), intent(in) :: arg
        real(_pr) :: res
        res = arg**2
    end function func

    function diff(arg, fptr) result(res)
        procedure(func), pointer :: fptr
        real(_pr), intent(in):: arg
        real(_pr) :: res
        res = (fptr(arg+del) - fptr(arg-del))/(2*del)
    end function diff

    function ddiff(arg, fptr) result(res)
        procedure(func), pointer :: fptr
        real(_pr), intent(in):: arg
        real(_pr) :: res
        res = (fptr(arg+del) - 2*fptr(arg) + fptr(arg-del))/(del**2)
    end function ddiff
end program draw

