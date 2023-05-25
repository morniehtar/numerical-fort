module global
    implicit none

    integer, parameter :: n = 100
    integer, parameter :: tries = 1000000 ! 1000000 is ok

    real(8), dimension(n) :: histo = 0

contains

    pure real(8) function gauss(arg, avg, sig)
        implicit none
        real(8), intent(in) :: arg, avg, sig
        real(8), parameter :: pi = 4*atan(1.d0)

        gauss = (dexp(-((arg-avg)**2)/(2.0d0*sig**2))/(sig*dsqrt(2.0d0*pi)))

    end function gauss

end module global


program main
    use global
    implicit none
    external chi
    integer :: i, j, k, vars, maxcal, ifail
    real(8) :: step, eps, fX(2), fZ(2), fMin

! Conducting experiment
    call random_seed()
    do i = 1, tries
        k = 0
        do j = 1, n
            if (flip()) k = k + 1
        end do
        histo(k) = histo(k) + 1
    end do

! Normalization
    histo = histo/tries

! Fitting
    vars = 2
    step = 0.05d0
    fX(1) = 100d0 ! expected average
    fX(2) = 25d0 ! dispersion
    eps = 1.d-12
    do while(.true.)
        maxcal = 100
        ifail = 0
        call SIMPLY(chi, vars, step, eps, fX, fZ, fMin, maxcal, ifail)
        !write(*,*) fZ, fMin, ifail
        if(ifail.eq.0) exit
        fX(1) = fZ(1)
        fX(2) = fZ(2)
    end do

! Output
    open(unit = 1, file = "plot.dat")
    do i = 1, n
        write(1,*) i, histo(i), gauss(dfloat(i), fZ(1), fZ(2))
    end do
    close(1)

    print*, "IFAIL status: ", ifail
    print*, "Mu: ", fZ(1)
    print*, "Sigma: ", fZ(2)
    print*, "Value: ", fMin

contains

    logical function flip()
        implicit none
        real :: r

        !call random_seed()
        call random_number(r)

        if (r .lt. 0.5d0) then
        flip = 0
        else
            flip = 1
        end if

    end function flip


end program main

subroutine chi(arg, ans)
    use global
    implicit none

    real(8), intent(in) :: arg(2)
    real(8), intent(out) :: ans
    integer :: j

    ans = 0.d0
    do j = 1, n
        ans = ans + (gauss(dfloat(j), arg(1), arg(2))-histo(j))**2
    end do

end subroutine chi
