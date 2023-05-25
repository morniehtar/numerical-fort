module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: fl = 0.d0, cl = 1.d0
    !Integration precision
    integer, parameter :: dots = 10
end module config

program output
    use config
    implicit none

    procedure(func), pointer :: fptr
    fptr => func

    call integ(fptr, fl, cl)


contains
    !y' = func(x, y)

    pure real(8) function func(x, y)
        implicit none
        real(8), intent(in) :: x, y
        func = 5*x**4
    end function func

    subroutine integ(fptr, st, ed)
        !Blue line
        use config
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: st, ed
        real(8) :: x0, x, y, k1, k2, k3, k4, step
        integer :: i

        step = (ed - st) / dots

        !Начальные условия
        x0 = 0.d0
        y = 0.d0

        open(unit = 1, file = "output.dat")
        do i = 1, dots+1
            x = x0 + (i-1)*step
            write(1,*) x, y, sltn(x)
            k1 = step*fptr(x, y)
            k2 = step*fptr(x+step*0.5d0, y+k1*step*0.5d0)
            k3 = step*fptr(x+step*0.5d0, y+k2*step*0.5d0)
            k4 = step*fptr(x+step, y+k3*step)
            y = y + (k1+2.d0*k2+2.d0*k3+k4)/6.d0
        end do
        close(1)

    end subroutine integ

    pure real(8) function sltn(x)
        !Red line
        implicit none
        real(8), intent(in) :: x
        sltn = x**5
    end function sltn

end program output
