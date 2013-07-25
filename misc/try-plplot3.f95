module plutil
contains

    subroutine gray_cmap(num_col)
        use plplot
        integer num_col
        real(kind=plflt) r(2), g(2), b(2), pos(2)
        logical rev(2)
        r(1) = 0.0
        g(1) = 0.0
        b(1) = 0.0
        r(2) = 1.0
        g(2) = 1.0
        b(2) = 1.0
        pos(1) = 0.0
        pos(2) = 1.0
        rev(1) = .false.
        rev(2) = .false.
        call plscmap1n(num_col)
        call plscmap1l(.true., pos, r, g, b, rev)
    end subroutine gray_cmap

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
        call gray_cmap(100)
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
