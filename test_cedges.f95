! gfortran -I /usr/lib/fortran/modules/plplot test-plplot4.f95 -lplplotf95d -lplplotd

program test_cedges
    use cedges
    real :: image(376,273)
    integer :: x,y
    logical :: flag

    image = 0
    image(100,100) = 1
    image(133,237) = 3
    x = 1
    y = 1
    flag = next_start(image,x,y)
    print *,x,y,flag
    if(x/=100 .or. y/=100) stop "assertion"
    image(x,y) = -1
    flag = next_start(image,x,y)
    print *,x,y,flag
    if(x/=133 .or. y/=237) stop "assertion"
    image(x,y) = -1
    flag = next_start(image,x,y)
    print *,x,y,flag

    image = 0
    do i=0,10
        image(100+i,100+i) = 1
    end do
    x = 1
    y = 1
    flag = next_start(image,x,y)
    do while(next_neighbor(image,x,y))
        print *,x,y
        image(x,y) = -1
    end do
    if (x/=110 .or. y/=110) stop "assertion"
end program test_cedges

