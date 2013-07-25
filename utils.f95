module utils
contains
    pure real function dist(u,v)
        real, intent(in) :: u(:),v(:)
        dist = sqrt(sum((u-v)**2))
    end function dist

    integer function randint(n)
        integer n
        real v
        call random_number(v)
        randint = modulo(floor(n*v),n)+1
    end function randint

    real function randunif()
        real value
        call random_number(value)
        randunif = value
    end function randunif

    real function logunif(lo,hi)
        real lo,hi
        logunif = lo+(hi-lo)*randunif()
    end function logunif

    real function logspace(i,n,lo,hi)
        integer i,n
        real lo,hi
        logspace = exp(((i-1) / real(n-1)) * (log(hi)-log(lo)) + log(lo))
    end function logspace

    real function randnormal()
        real x,y,s
        do
            x = 2*randunif()-1
            y = 2*randunif()-1
            s = x**2+y**2
            if (s<=1.0) exit
        end do
        randnormal = -sqrt(-log(s)/s)*x
    end function randnormal

    integer function imod(i,n)
        integer i,n
        imod = modulo(i-1,n)+1
    end function imod

    subroutine indent(n)
        do i=1,n
            write(*,fmt="(a)",advance="no") " "
        end do
    end subroutine indent

    subroutine perturb(v,s)
        real v(:),s
        integer i
        do i=1,size(v)
            v(i) = v(i) + s * randunif()
        end do
    end subroutine perturb

    subroutine hist(h,v)
        integer, allocatable :: h(:)
        integer :: v(:)
        integer i,c
        allocate(h(minval(v):maxval(v)))
        h = 0
        do i=1,size(v)
            c = v(i)
            h(c) = h(c) + 1
        end do
    end subroutine hist

    subroutine print_hist(v)
        integer :: v(:)
        integer, allocatable :: h(:)
        integer i,c,lo,hi
        lo = minval(v)
        hi = maxval(v)
        print *,"histogram",lo,hi
        allocate(h(lo:hi))
        h = 0
        do i=1,size(v)
            c = v(i)
            h(c) = h(c) + 1
        end do
        do i=minval(v),maxval(v)
            print *,"hist",i,h(i)
        end do
    end subroutine print_hist

    real function mindist(v,means)
        real, intent(in) :: v(:),means(:,:)
        real dists(size(means,1))
        integer i,k,d
        k = size(means,1)
        d = size(means,2)
        !$omp parallel shared(dists,means,d)
        !$omp do schedule(runtime)
        do i=1,k
            dists(i) = dist(v,means(i,:))
            if (dists(i)<0 .or. dists(i)>1e10) stop 'bad distance'
        end do
        !$omp end do
        !$omp end parallel
        mindist = minval(dists)
    end function mindist

    integer function argmindist(v,means)
        real, intent(in) :: v(:),means(:,:)
        real dists(size(means,1))
        integer i,k,d
        k = size(means,1)
        d = size(means,2)
        !$omp parallel shared(dists,means,d)
        !$omp do schedule(runtime)
        do i=1,k
            dists(i) = dist(v,means(i,:))
            if (dists(i)<0 .or. dists(i)>1e10) stop 'bad distance'
        end do
        !$omp end do
        !$omp end parallel
        argmindist = minloc(dists,1)
    end function argmindist

    subroutine argmindist2(i1,d1,i2,d2,v,means,k,d)
        integer, intent(out) :: i1,i2
        real, intent(out) :: d1,d2
        integer, intent(in) :: k,d
        real, intent(in) :: v(d),means(k,d)
        integer i
        real dists(k)

        !$omp parallel shared(dists,means,d)
        !$omp do schedule(runtime)
        do i=1,k
            dists(i) = dist(v,means(i,:))
            if (dists(i)<0 .or. dists(i)>1e10) stop 'bad distance'
        end do
        !$omp end do
        !$omp end parallel

        i1 = minloc(dists,1)
        d1 = dists(i1)
        dists(i1) = 1e30
        i2 = minloc(dists,1)
        d2 = dists(i2)
    end subroutine argmindist2


    real function fgetenv(name,default)
        character(*) name
        character(100) s
        real x,default
        integer status
        call get_environment_variable(name,s,status=status)
        if (status==0) then
            read(s,fmt=*) x
            fgetenv = x
        else
            fgetenv = default
        end if
        print *,name,"=",fgetenv
    end function fgetenv

    integer function igetenv(name,default)
        character(*) name
        character(100) s
        integer x,default
        integer status
        call get_environment_variable(name,s,status=status)
        if (status==0) then
            read(s,fmt=*) x
            igetenv = x
        else
            igetenv = default
        end if
        print *,name,"=",igetenv
    end function igetenv

    character(1000) function sgetenv(name,default)
        character(*) name,default
        character(1000) s
        integer status
        call get_environment_variable(name,s,status=status)
        if (status==0) then
            sgetenv = s
        else
            sgetenv = default
        end if
        print *,name,"=",sgetenv
    end function sgetenv

end module utils
