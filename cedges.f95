module cedges
    use reallocs
    use quicksort

    implicit none

    ! thresholds used by the Canny edge detector
    real :: frac = 0.3
    real :: tlow = 2.0
    real :: thigh = 4.0

    ! last noise estimate by the Canny edge detector
    real :: noise = -1

    ! max distance for polygonal approximation
    real :: maxdist = 1.5

    ! neighborhood table used for thinning
    integer :: thin_table(256) = &
        &[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, &
        &0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, &
        &0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        &0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, &
        &0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, &
        &0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, &
        &0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        &1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, &
        &0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, &
        &0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

    ! an enumeration of the x and y offsets of the neighbors of a pixel
    integer :: nx(8) = [1,1,0,-1,-1,-1,0,1]
    integer :: ny(8) = [0,1,1,1,0,-1,-1,-1]

contains
    subroutine gmask(mask,sigma)
        real, allocatable :: mask(:)
        real sigma,y
        integer range,i
        range = 1+floor(3.0*sigma)
        allocate(mask(2*range+1))
        do i=0,range
            y = exp(-i*i/2.0/sigma/sigma)
            mask(range+1+i) = y
            mask(range+1-i) = y
        end do
        mask = mask/sum(mask)
    end subroutine gmask

    subroutine convolve1d(signal,mask)
        real signal(:),mask(:),sigma,y
        real temp(size(signal))
        integer range,shift
        range = size(mask)/2
        temp = 0
        do shift=-range,range
            temp = temp + mask(range+1-shift) * eoshift(signal,shift,0.0,1)
        end do
        signal = temp
    end subroutine convolve1d

    subroutine gauss2d(image,sigma)
        real image(:,:)
        real sigma
        real, allocatable :: mask(:)
        integer i
        call gmask(mask,sigma)
        do i=1,size(image,1)
            call convolve1d(image(i,:),mask)
        end do
        do i=1,size(image,2)
            call convolve1d(image(:,i),mask)
        end do
    end subroutine gauss2d

    real function val(flag)
        logical flag
        if (flag) then
            val = 1.0
        else
            val = 0.0
        end if
    end function val

    subroutine nonmaxsup(out,gradm,gradx,grady)
        real out(:,:),gradm(:,:),gradx(:,:),grady(:,:)
        real dx,dy,ux,uy,vx,vy,c,u,d
        integer ax,ay,bx,by
        integer w,h,i,j
        w = size(gradm,1)
        h = size(gradm,2)
        out = 0
        do i=2,w-1
            do j=2,h-1
                dx = gradx(i,j)
                ux = abs(dx)
                dy = grady(i,j)
                uy = abs(dy)
                bx = sign(1.0,dx)
                by = sign(1.0,dy)
                ax = bx*val(ux>uy)
                ay = by*val(ux<=uy)
                if (ax/=0) then
                    vy = ux
                    vx = uy
                else 
                    vx = ux
                    vy = uy
                end if
                c = gradm(i,j)
                u = gradm(i-ax,j-ay)
                d = gradm(i-bx,j-by)
                if (vy*c<=(vx*d+(vy-vx)*u)) cycle
                u = gradm(i+ax,j+ay)
                d = gradm(i+bx,j+by)
                if (vy*c<(vx*d+(vy-vx)*u)) cycle
                out(i,j) = 1
            end do
        end do
    end subroutine nonmaxsup

    real function nonzero_fractile(gradm,frac)
        real gradm(:,:),frac
        real vector(count(gradm>0))
        vector = pack(gradm,gradm>0)
        call quicksortv(vector,size(vector))
        nonzero_fractile = vector(1+floor(frac*size(vector)))
    end function nonzero_fractile

    recursive subroutine masked_fill(out,gradm,tlow,x,y)
        real out(:,:),gradm(:,:),tlow
        integer x,y,i,j
        if (x<1.or.x>size(out,1).or.y<1.or.y>size(out,2)) return
        if (out(x,y)/=0) return
        if (gradm(x,y)<tlow) return
        out(x,y) = 1
        do i=-1,1
            do j=-1,1
                if(i==0.and.j==0) cycle
                call masked_fill(out,gradm,tlow,x+i,y+j)
            end do
        end do
    end subroutine masked_fill

    subroutine hysteresis_thresholding(out,gradm,tlow,thigh)
        real out(:,:),gradm(:,:),tlow,thigh
        integer i,j
        out = 0
        do i=1,size(gradm,1)
            do j=1,size(gradm,2)
                if (out(i,j)>0) cycle
                if (gradm(i,j)>thigh) call masked_fill(out,gradm,tlow,i,j)
            end do
        end do
    end subroutine hysteresis_thresholding

    subroutine thin(image)
        real image(:,:)
        integer b,dir,i,x,y
        integer, parameter :: OFF=0,ON=1,SKEL=2,DEL=3
        logical flag
        where (image>0)
            image = ON
        elsewhere
            image = OFF
        end where
        flag = .true.
        do while(flag)
            flag = .false.
            do dir=1,8
                do x=2,size(image,1)-1
                    do y=2,size(image,2)-1
                        if (image(x,y)/=ON) cycle
                        if (image(x+nx(dir),y+ny(dir))/=OFF) cycle
                        b = 0
                        do i=8,1,-1
                            b = 2*b
                            if (image(x+nx(i),y+ny(i))/=OFF) b=b+1
                        end do
                        if (thin_table(b+1)>0) then
                            image(x,y) = SKEL
                        else
                            image(x,y) = DEL
                            flag = .true.
                        end if
                    end do
                end do
                if (.not. flag) cycle
                where (image==DEL) image=OFF
            end do
        end do
        where (image==SKEL)
            image = 1
        elsewhere
            image = 0
        end where
    end subroutine thin

    subroutine canny(out,image)
        real image(:,:),out(:,:)
        real gradm(size(image,1),size(image,2))
        real gradx(size(image,1),size(image,2))
        real grady(size(image,1),size(image,2))
        integer i,j
        gradx = image - eoshift(image,-1,0.0,1)
        gradx(1,:) = 0
        grady = image - eoshift(image,-1,0.0,2)
        grady(:,1) = 0
        gradm = sqrt(gradx**2 + grady**2)
        call nonmaxsup(out,gradm,gradx,grady)
        call thin(out)
        where (out==0) gradm = 0
        noise = nonzero_fractile(gradm,frac)
        call hysteresis_thresholding(out,gradm,tlow*noise,thigh*noise)
    end subroutine canny

    integer function count_neighbors(image,x,y)
        real image(:,:)
        integer x,y,count,i
        count = 0
        do i=1,8
            if (image(x+nx(i),y+ny(i))>0) count = count+1
        end do
        count_neighbors = count
    end function count_neighbors

    logical function next_neighbor(image,x,y)
        real image(:,:)
        integer x,y,count,i
        count = 0
        do i=1,8
            if (image(x+nx(i),y+ny(i))>0) then
                x = x+nx(i)
                y = y+ny(i)
                next_neighbor = .true.
                return
            end if
        end do
        next_neighbor = .false.
    end function next_neighbor

    logical function next_start(image,x,y)
        real image(:,:)
        integer x,y
        do while(x<=size(image,1))
            do while(y<=size(image,2))
                if (image(x,y)>0) then
                    next_start = .true.
                    return
                end if
                y = y+1
            end do
            x = x+1
            y = 1
        end do
        next_start = .false.
    end function next_start

    subroutine reverse_rows(a)
        integer a(:,:)
        integer temp(size(a,1),size(a,2))
        temp(size(a,1):1:-1,:) = a(:,:)
        a = temp
    end subroutine reverse_rows

    ! Trace the lines in a thinned image and store the traced
    ! coordinates in the points array.
    ! Breaks between traces are not explicitly indicated; they
    ! are simply wherever there is more than 1 pixel difference
    ! between one pixel and the next.

    subroutine chains(points,image)
        integer, allocatable :: points(:,:)
        integer npoints
        real image(:,:)
        integer last,x,y,xs,ys
        image(1,:) = 0
        image(:,1) = 0
        image(size(image,1),:) = 0
        image(:,size(image,2)) = 0
        allocate(points(1000000,2))
        npoints = 0
        xs = 1
        ys = 1
        do while(next_start(image,xs,ys))
            x = xs
            y = ys
            image(x,y) = -1
            print *,"x,y",x,y
            npoints = npoints+1
            if (npoints>size(points,1)) return
            points(npoints,:) = [x,y]
            last = npoints
            do while(next_neighbor(image,x,y))
                npoints = npoints+1
                if (npoints>size(points,1)) return
                points(npoints,:) = [x,y]
                image(x,y) = -1
            end do
            ! call reverse_rows(points(last:npoints,:))
            do while(next_neighbor(image,points(npoints,1),points(npoints,2)))
                npoints = npoints+1
                if (npoints>size(points,1)) return
                points(npoints,:) = [x,y]
                image(x,y) = -1
            end do
        end do
        call realloc(points,npoints,2)
    end subroutine chains

    real function point_line_dist(p,a,b)
        integer p(2),a(2),b(2)
        real delta(2),normal(2)
        real mag
        delta = b-a
        mag = norm(delta)
        if (mag<1e-4) then
            point_line_dist = inorm(a-p)
        else
            normal = [-delta(2),delta(1)]/mag
            point_line_dist = abs(dot_product(normal,p) - dot_product(normal,a))
        end if
    end function point_line_dist

    real function norm(x)
        real x(:)
        norm = sqrt(sum(x**2))
    end function norm

    real function inorm(x)
        integer x(:)
        inorm = sqrt(1.0*sum(x**2))
    end function inorm

    recursive subroutine poly_approx0(polys,npolys,chain)
        integer polys(:,:)
        integer npolys
        integer chain(:,:)
        real dists(size(chain,1))
        integer split,i
        if (size(dists)<2) return
        do i=1,size(dists)
            dists(i) = point_line_dist(chain(i,:),chain(1,:),chain(size(chain,1),:))
        end do
        split = maxloc(dists,1)
        if (dists(split)<maxdist) then
            npolys = npolys + 1
            polys(npolys,1:2) = chain(1,:)
            polys(npolys,3:4) = chain(size(chain,1),:)
            return
        end if
        if (split>1) call poly_approx0(polys,npolys,chain(1:split,:))
        if (split<size(chain,1)) call poly_approx0(polys,npolys,chain(split:size(chain,1),:))
    end subroutine poly_approx0

    ! Approximate the pixel chains in the chain array with polygons.
    ! Deviation from the true chain may not be larger than maxdist.

    subroutine poly_approx(polys,chain)
        integer, allocatable :: polys(:,:)
        integer npolys
        integer chain(:,:)
        integer start,last,here,i
        real d
        npolys = 0
        allocate(polys(1000000,4))
        polys = 0
        last = 1
        do i=2,size(chain,1)
            if (i==size(chain) .or. abs(chain(i,1)-chain(i-1,1))>2 .or. &
                &abs(chain(i,2)-chain(i-1,2))>2) then
                ! npolys = npolys + 1
                ! polys(npolys,:) = [chain(last,1),chain(last,2),chain(i-1,1),chain(i-1,2)]
                call poly_approx0(polys,npolys,chain(last:i-1,:))
                last = i
            end if
        end do
        call realloc(polys,npolys,4)
    end subroutine poly_approx

    subroutine cedges_test
        use jpegio
        real, allocatable :: image(:,:), raw(:,:)
        real, allocatable :: out(:,:)
        integer, allocatable :: points(:,:), polys(:,:)
        integer npoints,npolys,i
        call read_jpeg2("test.jpg",raw)
        raw = raw/maxval(raw)
        allocate(image(size(raw,1),size(raw,2)))
        allocate(out(size(image,1),size(image,2)))
        do i=1,10
            image = raw
            call gauss2d(image,3.0)
            call canny(out,image)
            ! out = 255*out/maxval(out)
            ! call write_jpeg2("out.jpg",out)
            call chains(points,out)
            call poly_approx(polys,points)
        end do
        do i=1,size(polys,1)
            print *,i,polys(i,:)
        end do
    end subroutine cedges_test

    subroutine get_points(points,fname)
        use jpegio
        character(*) fname
        integer, allocatable :: points(:,:)
        real, allocatable :: image(:,:),edges(:,:)
        integer npoints
        call read_jpeg2(fname,image)
        image = image / maxval(image)
        image(:,:) = image(:,size(image,2):1:-1)
        call gauss2d(image,3.0)
        allocate(edges(size(image,1),size(image,2)))
        call canny(edges,image)
        call chains(points,edges)
    end subroutine get_points

    subroutine get_segments(polys,fname)
        use jpegio
        character(*) fname
        integer, allocatable :: polys(:,:),points(:,:)
        real, allocatable :: image(:,:),edges(:,:)
        integer npoints,npolys
        call read_jpeg2(fname,image)
        image = image / maxval(image)
        image(:,:) = image(:,size(image,2):1:-1)
        allocate(edges(size(image,1),size(image,2)))
        call gauss2d(image,3.0)
        call canny(edges,image)
        call chains(points,edges)
        call poly_approx(polys,points)
    end subroutine get_segments
end module cedges

