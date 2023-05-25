module config
    implicit none

! Potential type (1 for linear, 2 for quadratic)
    integer :: pt = 2
! Potential scaling factor
    real, parameter :: alpha=1e6
! Plot smoothness
    integer, parameter :: dots=1e3
! Eigenvalues error
    real, parameter :: eps=1e-6
! Voevdn iterations number
    integer, parameter :: itr=10000
! Dimensions
    integer, parameter :: ndim = 16

    real, parameter :: pi = 4*atan(1e0)
end module config

program main
    use config
    implicit none

    real, dimension(ndim, ndim):: V, eigen
    integer :: i, j
    real :: tausq

! Perturbation matrix, n by n
    do i=1, ndim
        do j=1, ndim
            if(i==j) then
                if (pt == 1) then
                    V(i,j)=alpha/2e0+(pi**2*i**2)/2e0
                elseif (pt == 2) then
                    V(i,j)=alpha*(2*i*i*pi*pi-3.0d0)/(12.0d0*i*i*pi*pi)+pi*pi*i*i/2.0d0
                end if
            else
                if (pt == 1) then
                    V(i,j)=(alpha/(pi**2))*(((-1e0)**(i-j)-1e0)/((i-j)**2)+(1e0-(-1e0)**(i+j))/((i+j)**2))
                elseif (pt == 2) then
                    V(i,j)=(alpha/(pi*pi))*(((-1.0d0)**(i-j))/((i-j)*(i-j))-((-1.0d0)**(i+j))/((i+j)*(i+j)))
                end if
            end if
        end do
    end do

    eigen = 0d0
    tausq = 0d0
    call voevdn(V, eigen, ndim, itr, eps, 0, tausq)

    write (*, "(4a16)") "Unperturbed:", "Perturbed:", "Difference:", "Relative:"
    do i = 1, ndim
        write(*, "(i2.2, 1000f16.5)") i, 0.5d0*pi**2*i**2, V(i, i), V(i, i) - 0.5d0*pi**2*i**2, &
        (V(i, i) - 0.5d0*pi**2*i**2)/(0.5d0*pi**2*i**2)
    end do

    !call drawMtx(eigen)
    !print*, "Off-diagonal sum: ", tausq

    do i = 1, ndim
        call drawWF(eigen)
    end do

contains

    subroutine drawMtx(arg)
        real, dimension(ndim, ndim), intent(in) :: arg
        integer :: i_dr, j_dr

        do i_dr=1, ndim
            write(*, "(1000e13.5)") ( arg(i_dr, j_dr), j_dr=1 , ndim )
        end do

    end subroutine

    subroutine drawWF(c)
        implicit none
        real, dimension(ndim, ndim), intent(in) :: c
        real :: x, psi(ndim)
        character(2) :: chDim, chAlph
        integer :: i_wf, j_wf, lv

    ! File naming
        write(chDim, '(I2.2)') ndim
        i_wf = alpha
        j_wf = 0
        do while (i_wf > 1)
            i_wf = i_wf / 10
            j_wf = j_wf + 1
        end do
        write(chAlph, '(I2.2)') j_wf
        if (pt == 1) then
            open(1, file="./wf01data/wf_n"//trim(chDim)//"_a10e"//trim(chAlph)//".dat")
        elseif (pt == 2) then
            open(1, file="./wf02data/wf_n"//trim(chDim)//"_a10e"//trim(chAlph)//".dat")
        end if

    ! WF drawing

        do i_wf = 0, dots
            x = 0e0 + 1e0*i_wf/dots
            psi=0e0
            do lv = 1, ndim
                do j_wf = 1, ndim
                    psi(lv)=psi(lv)+c(j_wf, lv)*sin(pi*j_wf*x)*sqrt(2e0)
                end do
            end do
            write(1,*) x, (psi(lv), lv = 1, ndim)
        end do
        close(1)
    end subroutine drawWF

end program main




