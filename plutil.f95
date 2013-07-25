module plutil
    use plplot
    implicit none
contains

    real function distance(u,v)
        real u(:),v(:)
        distance = sqrt(sum((u-v)**2))
    end function distance

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

    subroutine pluimage(image)
        use plplot
        real :: image(:,:)
        real(plflt), allocatable :: zs(:,:)
        real(plflt) :: clevel(100)
        real(plflt) :: w,h
        integer i
        allocate(zs(size(image,1),size(image,2)))
        ! zs(:,:) = image(:,size(image,2):1:-1)
        zs = image
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

    subroutine plusegments(segs)
        integer segs(:,:),i
        real(plflt) xs(2),ys(2)
        if (size(segs,2)/=4) stop "bad argument to plusements"
        do i=1,size(segs,1)
            xs(:) = [segs(i,1),segs(i,3)]
            ys(:) = [segs(i,2),segs(i,4)]
            call plline(xs,ys)
        end do
    end subroutine plusegments

    subroutine pluchain(chain)
        integer chain(:,:),i,last
        real(plflt) xs(100000),ys(100000)
        if (size(chain,2)/=2) stop "bad argument to pluchain"
        last = 1
        do i=2,size(chain,1)
            if (i==size(chain) .or. abs(chain(i,1)-chain(i-1,1))>2 .or. &
                &abs(chain(i,2)-chain(i-1,2))>2) then
                print *,i,last,size(xs,1)
                xs(:i-last) = chain(last:i-1,1)
                ys(:i-last) = chain(last:i-1,2)
                call plline(xs(:i-last),ys(:i-last))
                last = i
            end if
        end do
    end subroutine pluchain
end module plutil

