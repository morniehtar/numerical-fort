module config
    implicit none

! Try seven parameters
    integer, parameter :: N = 4
    real(8), parameter :: g = 2d0

    real(8), parameter :: bt = 0d0
    real(8), parameter :: tp = 7d0
    integer, parameter :: dots = 1d3

    integer, parameter :: msteps = 1000
    real(8), parameter :: eps = 1d-12
    real(8), parameter :: stp = 5d-2

end module config



program main
    use config
    implicit none
    external func
    external fct

    integer :: ifail, i, j, m
    real(8), dimension(N) :: x, z
    real(8) :: fe, fct, y, r, norm


! --------------------
    do i = 1, N
        x(i)=3d0
    enddo

    do while(.true.)
        ifail=0
        call simply(func, N, stp, eps, x, z, fe, msteps, ifail)
        if (ifail.eq.0) exit

        do i = 1, N
            x(i) = z(i)
        enddo
    enddo

    print*, "Parameters: "
    do i = 1, N
        print*, z(i)
    enddo
    print*, "Minimal value: ", fe

    open(2, file='plot.dat')
    do i = 0, dots
        r=0d0
        norm=0d0
        r = bt + (tp-bt)*i/dots
        do j=1,N-1
            y = y+z(j)*(r**j)
        enddo
        y = y*dexp(-r*z(N)/2d0)
        do j = 1, N-1
            do m = 1, N-1
                norm = norm+z(j)*z(m)*fct(m+j)/(x(N)**(j+m+1))
            enddo
        enddo
        y = y/sqrt(norm)
        write(2,*) r, y
    enddo
    close(2)


end program main

subroutine func(x, res)
    use config
    implicit none
    real(8), dimension(N), intent(in) :: x
    real(8), intent(out) :: res
    real(8) :: denom
    integer :: i_fu, j_fu

    real(8) :: fct
    external fct

    res=0.0d0
    denom=0.0d0
    do i_fu = 1, N-1
        do j_fu = 1, N-1
            res=res+x(j_fu)*x(i_fu)*(0.5d0*fct(i_fu+j_fu-2)/x(N)**(i_fu+j_fu-1)*(i_fu*j_fu-0.5*(i_fu+j_fu)*(i_fu+j_fu-1)&
                +0.25d0*(i_fu+j_fu)*(i_fu+j_fu-1))-g*fct(i_fu+j_fu-1)/(x(N)+1)**(i_fu+j_fu))
            denom=denom+x(i_fu)*x(j_fu)*fct(j_fu+i_fu)/(x(N)**(i_fu+j_fu+1))
        enddo
    enddo
    res=res/denom
end subroutine

pure real(8) function fct(num)
    integer, intent(in) :: num
    integer :: i_fct

    fct=1d0
    if (num.le.0) return

    i_fct=num
    do i_fct = num, 1, -1
        fct = fct*i_fct
    end do
end function fct
