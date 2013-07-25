module reallocs
    interface realloc
        module procedure realloc1i, realloc1r, realloc2i, realloc2r,&
            &realloc1l,realloc2l,realloc1i2,realloc2i2
    end interface
contains
    subroutine realloc1i(a,n)
        integer, allocatable :: a(:),temp(:)
        integer n,n0
        allocate(temp(n))
        n0 = min(size(a),n)
        temp(1:n0) = a(1:n0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc1i
    subroutine realloc2i(a,n,m)
        integer, allocatable :: a(:,:),temp(:,:)
        integer n,m,n0,m0
        allocate(temp(n,m))
        n0 = min(size(a,1),n)
        m0 = min(size(a,2),m)
        temp(1:n0,1:m0) = a(1:n0,1:m0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc2i

    subroutine realloc1i2(a,n)
        integer(2), allocatable :: a(:),temp(:)
        integer n,n0
        allocate(temp(n))
        n0 = min(size(a),n)
        temp(1:n0) = a(1:n0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc1i2
    subroutine realloc2i2(a,n,m)
        integer(2), allocatable :: a(:,:),temp(:,:)
        integer n,m,n0,m0
        allocate(temp(n,m))
        n0 = min(size(a,1),n)
        m0 = min(size(a,2),m)
        temp(1:n0,1:m0) = a(1:n0,1:m0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc2i2

    subroutine realloc1l(a,n)
        logical, allocatable :: a(:),temp(:)
        integer n,n0
        allocate(temp(n))
        n0 = min(size(a),n)
        temp(1:n0) = a(1:n0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc1l
    subroutine realloc2l(a,n,m)
        logical, allocatable :: a(:,:),temp(:,:)
        integer n,m,n0,m0
        allocate(temp(n,m))
        n0 = min(size(a,1),n)
        m0 = min(size(a,2),m)
        temp(1:n0,1:m0) = a(1:n0,1:m0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc2l

    subroutine realloc1r(a,n)
        real, allocatable :: a(:),temp(:)
        integer n,n0
        allocate(temp(n))
        n0 = min(size(a),n)
        temp(1:n0) = a(1:n0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc1r
    subroutine realloc2r(a,n,m)
        real, allocatable :: a(:,:),temp(:,:)
        integer n,m,n0,m0
        allocate(temp(n,m))
        n0 = min(size(a,1),n)
        m0 = min(size(a,2),m)
        temp(1:n0,1:m0) = a(1:n0,1:m0)
        call move_alloc(from=temp,to=a)
    end subroutine realloc2r
end module reallocs
