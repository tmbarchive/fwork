module plutil
contains
    subroutine pluimage(zs)
        use plplot
        real(plflt) :: zs(:,:)
        real(plflt) :: clevel(100)
        real(plflt) :: w,h
        integer i
        w = size(zs,1)
        h = size(zs,2)
        forall (i=1:100) clevel(i) = (i-1)/99.0
        call plenv(1.0_plflt,w,1.0_plflt,h,1,1)
        call plshades(z=zs,defined="all",&
            &xmin=1.0_plflt,xmax=w,&
            &ymin=1.0_plflt,ymax=h,&
            &clevel=clevel,fill_width=0,cont_color=0,cont_width=0)
    end subroutine pluimage
end module plutil

program testpl
    use plplot
    use plutil
    real(plflt) :: zs(100,100)
    do i=1,100
        do j=1,100
            zs(i,j) = i*j/10000.0
        end do
    end do
    call plsdev('gcw')
    call plinit
    call pluimage(zs)
    call plend
end program testpl
