module search
    implicit none
    private

    abstract interface
        real(8) function func(x)
        real(8), intent(in) :: x
        end function func
    end interface

    abstract interface
        subroutine srch(fptr, bot, top, arr)
        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: bot, top
        real(8), dimension(:), intent(inout), allocatable :: arr
        end subroutine srch
    end interface

! Root search method
    procedure(srch), pointer :: getRoot => secant

! Root search limiter (affects number of correct digits in stationary levels)
    real(8), parameter :: eps = 1d-9 ! 1d-16 is machine epsilon

    public :: getRoot
contains

    subroutine bisection(fptr, bot, top, arr)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: bot, top
        real(8), dimension(:), intent(inout), allocatable :: arr
        real(8) :: md, bt, tp
        integer :: j

        if (fptr(bot)*fptr(top).gt.0) return
        call ext(arr)


        bt = bot
        tp = top
        md = (top + bot) / 2
        j = 0
        do while ((abs(fptr(md)).gt.eps).and.((tp-bt).gt.eps))
            md = (tp + bt) / 2
            if (fptr(bt) * fptr(md) .gt. 0) then
                bt = md
            else if (fptr(bt) * fptr(md) .lt. 0) then
                tp = md
            end if
            j = j + 1
        end do
        arr(size(arr)) = md

        print *, "(k =", size(arr), ") E =", md, j , "steps"
    end subroutine bisection

    subroutine secant(fptr, bot, top, arr)
        ! Modified to exclude dupes
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: bot, top
        real(8), dimension(:), intent(inout), allocatable :: arr
        real(8) :: bt, md, tp
        integer :: j, i

        if (fptr(bot)*fptr(top).gt.0) return
        call ext(arr)

        bt = bot
        tp = top
        md = tp - fptr(tp)*(tp-bt)/(fptr(tp)-fptr(bt))
        j = 0
        do while (abs(fptr(md)).gt.eps)
            md = tp - fptr(tp)*(tp-bt)/(fptr(tp)-fptr(bt))
            bt = tp
            tp = md
            j= j + 1
        end do

        do i = 1, size(arr)
            if (((abs(arr(i)-md)).lt.eps).or.(md.ge.0)) then
                call cut(arr)
                return
            end if
        end do
        arr(size(arr)) = md

        print *, "(k =", size(arr), ") E =", md, j , "steps"
    end subroutine secant

    subroutine ext(arr)
        implicit none
        real(8), dimension(:), intent(inout), allocatable :: arr

        real(8), dimension(size(arr)) :: tmp
        integer :: j, k

        k = size(arr)
        do j = 1, k
            tmp(j) = arr(j)
        end do

        deallocate(arr)
        allocate(arr(k+1))

        do j = 1, k
            arr(j) = tmp(j)
        end do
    end subroutine ext

    subroutine cut(arr)
        implicit none
        real(8), dimension(:), intent(inout), allocatable :: arr

        real(8), dimension(size(arr)) :: tmp
        integer :: j, k

        k = size(arr)
        do j = 1, k
            tmp(j) = arr(j)
        end do

        deallocate(arr)
        allocate(arr(k-1))

        do j = 1, k-1
            arr(j) = tmp(j)
        end do
    end subroutine cut

end module search
