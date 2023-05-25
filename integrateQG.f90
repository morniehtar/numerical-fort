module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: min = 0, max = 1
    !Integration step
    real(8), parameter :: step = 1d-6
    !Threads number
    integer, parameter :: thr = 8
end module config


program hello
    use config
    use omp_lib
    implicit none

    print *, "QG6: ", QG6(func, min, max, step)
    print *, "QG8: ", QG8(func, min, max, step)
    print *, "QG10:", QG10(func, min, max, step)

contains
    pure real(8) function func(x)
        real(8), intent(in) :: x
        func = x**2
    end function func

    real(8) function QG6(func, min, max, step)

        external func
        real(8) :: func
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(6) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.467913934572691d0, 0.467913934572691d0, 0.360761573048139d0, &
                0.360761573048139d0, 0.171324492379170d0, 0.171324492379170d0/
        data a /0.238619186083197d0,-0.238619186083197d0, 0.661209386466265d0, &
                -0.661209386466265d0, 0.932469514203152d0, -0.932469514203152d0/

        n = (max-min)/step
        half = step/2
        QG6=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG6) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG6)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 6
                    next = half*w(j)*func(a(j)*half+mid) - er
                    tmp = QG6 + next
                    er = (tmp - QG6) - next
                    QG6 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG6

    real(8) function QG8(func, min, max, step)

        external func
        real(8) :: func
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(8) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.362683783378362d0, 0.362683783378362d0, 0.313706645877887d0, 0.313706645877887d0, &
                0.222381034453374d0, 0.222381034453374d0, 0.101228536290376d0, 0.101228536290376d0/
        data a /0.183434642495650d0, -0.183434642495650d0, 0.525532409916329d0, -0.525532409916329d0, &
                0.796666477413627d0, -0.796666477413627d0, 0.960289856497536d0, -0.960289856497536d0/

        n = (max-min)/step
        half = step/2
        QG8=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG8) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG8)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 8
                    next = half*w(j)*func(a(j)*half+mid) - er
                    tmp = QG8 + next
                    er = (tmp - QG8) - next
                    QG8 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG8

    real(8) function QG10(func, min, max, step)

        external func
        real(8) :: func
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(10) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.066671344308688d0, 0.149451349150581d0, 0.219086362515982d0, &
                0.269266719309996d0, 0.295524224714753d0, 0.295524224714753d0, &
                0.269266719309996d0, 0.219086362515982d0, 0.149451349150581d0, &
                0.066671344308688d0/
        data a /0.97390652851717174d0, 0.86506336668898454d0, 0.67940956829902444d0, &
                0.43339539412924721d0, 0.14887433898163122d0, -0.14887433898163119d0, &
                -0.43339539412924721d0, -0.67940956829902444d0, -0.86506336668898454d0, &
                -0.97390652851717174d0/

        n = (max-min)/step
        half = step/2
        QG10=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG10) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG10)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 10
                    next = half*w(j)*func(a(j)*half+mid) - er
                    tmp = QG10 + next
                    er = (tmp - QG10) - next
                    QG10 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG10

end program
