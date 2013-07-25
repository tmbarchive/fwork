! gfortran -I /usr/lib/fortran/modules/plplot test-plplot4.f95 -lplplotf95d -lplplotd

program testpl
    use plplot
    use plutil
    use cedges
    use jpegio

    real, allocatable :: image(:,:), raw(:,:)
    real, allocatable :: out(:,:)
    integer, allocatable :: points(:,:), polys(:,:)
    integer npoints,i,npolys
    integer last

    call read_jpeg2("test.jpg",image)
    image = image / maxval(image)
    image(:,:) = image(:,size(image,2):1:-1)
    allocate(out(size(image,1),size(image,2)))
    call gauss2d(image,3.0)
    call canny(out,image)
    call chains(points,out)
    call poly_approx(polys,points)
    where (out/=0) out = 1

    call plsdev('gcw')
    call plinit
    call plenv(1.0_plflt,1.0_plflt*size(image,1),1.0_plflt,1.0_plflt*size(image,2),1,0)
    ! call pluimage(out)
    print *,"npoints",npoints
    print *,"npolys",npolys
    call plcol0(3)
    call plusegments(polys(:npolys,:))
    call plend
end program testpl
