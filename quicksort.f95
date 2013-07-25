module quicksort

contains

    subroutine swap(a,b)
        integer, intent(inout) :: a,b
        integer temp
        temp = a
        a = b
        b = temp
    end subroutine swap

    recursive subroutine quicksort_r(list,values,start,end,n)
        integer, intent(inout), dimension(0:n-1) :: list
        real, intent(in), dimension(0:n-1) :: values
        integer, intent(in) :: start, end, n
        integer :: lo, hi, split1, split2
        if (start>=end-1) return
        pivot = values(list((start+end-1)/2))
        lo = start
        hi = end
        do
            do while (lo<hi)
                if (.not. values(list(lo))<pivot) exit
                lo = lo + 1
            end do
            do while (lo<hi)
                if (.not. values(list(hi-1))>=pivot) exit
                hi = hi - 1
            end do
            if (lo>=hi-1) exit
            call swap(list(lo),list(hi-1))
            lo = lo + 1
            hi = hi - 1
        end do
        split1 = lo
        hi = end
        do
            do while (lo<hi)
                if (.not. values(list(lo))==pivot) exit
                lo = lo + 1
            end do
            do while (lo<hi)
                if (.not. values(list(hi-1))>pivot) exit
                hi = hi - 1
            end do
            if (lo>=hi-1) exit
            call swap(list(lo),list(hi-1))
            lo = lo + 1
            hi = hi - 1
        end do
        split2 = lo
        call quicksort_r(list,values,start,split1,n)
        call quicksort_r(list,values,split2,end,n)
    end subroutine quicksort_r

    subroutine quicksort0(list,values,n)
        integer, intent(out), dimension(0:n-1) :: list
        real, intent(in), dimension(0:n-1) :: values
        integer, intent(in) :: n
        forall (i=0:n-1) list(i) = i
        call quicksort_r(list,values,0,n,n)
    end subroutine quicksort0

    subroutine quicksort1(list,values,n)
        integer, intent(out), dimension(n) :: list
        real, intent(in), dimension(n) :: values
        integer, intent(in) :: n
        forall (i=1:n) list(i) = i-1
        call quicksort_r(list,values,0,n,n)
        forall (i=1:n) list(i) = list(i)+1
    end subroutine quicksort1

    subroutine quicksortv(values,n)
        integer list(n)
        real values(n),temp(n)
        call quicksort1(list,values,n)
        temp(:) = values(list(:))
        values(:) = temp(:)
    end subroutine quicksortv

    subroutine permutation1(list,n)
        integer n,list(n)
        real temp(n)
        call random_number(temp)
        call quicksort1(list,temp,n)
    end subroutine permutation1

    subroutine shuffle(values,n)
        integer n
        real values(n),temp(n)
        integer list(n)
        call random_number(temp)
        call quicksort1(list,temp,n)
        temp(:) = values(list(:))
        values(:) = temp(:)
    end subroutine shuffle

    subroutine shufflei(values,n)
        integer n
        real temp(n)
        integer values(n)
        integer list(n)
        integer tempi(n)
        call random_number(temp)
        call quicksort1(list,temp,n)
        tempi(:) = values(list(:))
        values(:) = tempi(:)
    end subroutine shufflei

    real function fractile_modify(a,f)
        real a(:)
        call quicksortv(a,size(a))
        index = floor(f*size(a)+1)
        fractile_modify = a(index)
    end function fractile_modify

    subroutine test
        real, dimension(5) :: a
        integer, dimension(5) :: l
        forall (i=1:5) a(i) = floor(10*sin(i*3.0))
        print *,a(:)
        call quicksort1(l,a,5)
        print *,a(l(:))
        print *,l
    end subroutine test
end module quicksort
