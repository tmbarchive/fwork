module rast
    use reallocs
    implicit none

    type segment
        real x0,y0,x1,y1
        real a
    end type segment

    type point
        real x,y
        real a
    end type point

    type matchlist
        integer(2), allocatable :: pairs(:,:)
        real, allocatable :: weights(:)
    end type matchlist

    real, parameter :: pi = 3.14159265358979323846

    integer :: verbosity = 10
contains
    elemental function hypot(a,b)
        real hypot
        real, intent(in) :: a,b
        hypot = sqrt(a**2+b**2)
    end function hypot

    elemental function heaviside(a)
        real heaviside
        real, intent(in) :: a
        if (a<0) then
            heaviside = 0.0
        else
            heaviside = a
        end if
    end function heaviside

    pure real function delta(iv)
        real, intent(in) :: iv(2)
        delta = iv(2)-iv(1)
    end function delta

    pure function deltav(iv)
        real, intent(in) :: iv(:,:)
        real deltav(size(iv,1))
        integer i
        forall (i=1:size(iv,1)) deltav(i) = iv(i,2)-iv(i,1)
    end function deltav

    pure real function mid(iv)
        real, intent(in) :: iv(2)
        mid = (iv(1)+iv(2))/2
    end function mid

    real function segment_point_dist(x0,y0,x1,y1,x,y)
        real x,y,x0,y0,x1,y1
        real dx,dy,nx,ny,r,l
        dx = x1-x0
        dy = y1-y0
        r = hypot(dx,dy)
        dx = dx / r
        dy = dy / r
        l = dx*x+dy*y
        if (l<dx*x0+dy*y0) then
            segment_point_dist = hypot(x-x0,y-y0)
        else if(l>dx*x1+dy*y1) then
            segment_point_dist = hypot(x-x1,y-y1)
        else
            segment_point_dist = abs((-dy*x+dx*y) - (-dy*x0+dx*y0))
        end if
    end function segment_point_dist

    real function score(d,eps)
        real d,eps
        if (d>eps) then
            score = 0
        else
            score = 1.0 - (d*d)/(eps*eps)
        end if
    end function score

    function d_score(d,r,eps)
        real, dimension(2) :: d_score
        real d,r,eps
        real d_lo,d_hi
        d_lo = heaviside(d-r)
        d_hi = d+r
        d_score(1) = score(d_hi,eps)
        d_score(2) = score(d_lo,eps)
    end function d_score

    function a_score(da,ra,aeps)
        real, dimension(2) :: a_score
        real da,ra,aeps,a,a_lo,a_hi
        if (ra>pi) then
            a_score = [0.0,1.0]
        else
            a = da
            do while (a>pi/2)
                a = a - pi
            end do
            do while (a<-pi/2)
                a = a + pi
            end do
            a = abs(a)
            a_lo = heaviside(a-ra)
            a_hi = heaviside(a+ra)
            a_score(1) = score(a_hi,aeps)
            a_score(2) = score(a_lo,aeps)
        end if
    end function a_score

    subroutine split(subtr,tr)
        real, dimension(4,2) :: tr
        real, allocatable, dimension(:,:,:) :: subtr
        real deltas(4)
        integer i,index
        forall (i=1:4) deltas(i) = delta(tr(i,:))
        deltas(3:4) = deltas(3:4) * 200.0
        if (minval(deltas)<0) stop "error"
        index = maxloc(deltas,1)
        allocate(subtr(2,4,2))
        subtr(1,:,:) = tr
        subtr(2,:,:) = tr
        subtr(1,index,2) = mid(tr(index,:))
        subtr(2,index,1) = mid(tr(index,:))
    end subroutine split

    subroutine decimate(eml,ml)
        type(matchlist) eml,ml
        integer i,j,n
        n = count(ml%pairs(:,1)>0)
        allocate(eml%pairs(n,2))
        allocate(eml%weights(n))
        eml%pairs(:,1) = pack(ml%pairs(:,1),ml%pairs(:,1)>0)
        eml%pairs(:,2) = pack(ml%pairs(:,2),ml%pairs(:,1)>0)
        eml%weights(:) = pack(ml%weights(:),ml%pairs(:,1)>0)
    end subroutine decimate

    ! Evaluate points against points.  The model consists of points,
    ! the image consists of segments.

    function evaluate_pp(points,ipoints,ml,tr,eps,aeps) result(score)
        type(point) points(:)
        type(point) ipoints(:)
        type(matchlist) ml
        real tr(4,2),eps,aeps

        type(segment) seg
        type(point) pt,ipt
        integer last_point,i,j
        real rfixed,rscaled,rtotal,tx,ty,tangle,tcos,tsin,tscale,px,py,pa
        real d,d_scores(2),a_scores(2),total(2),score(2)

        ! the actual transformations
        tx = mid(tr(1,:))
        ty = mid(tr(2,:))
        tangle = mid(tr(3,:))
        tscale = mid(tr(4,:))
        tcos = tscale * cos(tangle)
        tsin = tscale * sin(tangle)

        ! uncertainty that's absolute (independent of position)
        rfixed = hypot(delta(tr(1,:)),delta(tr(2,:)))

        ! uncertainty that scales with the distance of the vector from the origin
        rscaled = hypot(delta(tr(3,:)),delta(tr(4,:)))

        last_point = -1
        total = 0
        do i=1,size(ml%pairs,1)
            j = ml%pairs(i,1)
            if (j==-1) cycle
            if (last_point/=j) then
                pt = points(ml%pairs(i,1))
                px = tcos * pt%x - tsin * pt%y + tx
                py = tsin * pt%x + tcos * pt%y + ty
                pa = pt%a + tangle
                last_point = j
                rtotal = rfixed + hypot(px,py) * rscaled
            end if
            ipt = ipoints(ml%pairs(i,2))
            d = hypot(pt%x-ipt%x,pt%y-ipt%y)
            d_scores = d_score(d,rtotal,eps) * ml%weights(i)
            if (d_scores(2)<1e-6) then
                ml%pairs(i,1) = -1
            else
                total = total + d_scores
            end if
        end do
        score = total
    end function evaluate_pp

    subroutine search_pp(points,ipoints,ml,tr,eps,aeps,minq)
        real, parameter :: mindelta(4) = [1.0,1.0,0.01,0.002]
        type(point) points(:),ipoints(:),ipt,pt
        type(matchlist) ml
        real tr(4,2),eps,aeps,minq
        real best_score(2),best_tr(4,2)
        integer counts

        counts = 0
        print *,"start",size(ml%pairs,1),size(ml%pairs,2),&
            &minval(ml%weights),maxval(ml%weights)
        best_score = [0.0,0.0]
        call recursive(ml,tr,1)
    contains
        recursive subroutine recursive(ml,tr,depth)
            integer depth
            real tr(4,2)
            real score(2)
            real, allocatable :: subtr(:,:,:)
            type(matchlist) ml,eml
            integer i

            if (depth>100) then
                print *,"maxdepth exceeded"
                return 
            end if

            score = evaluate_pp(points,ipoints,ml,tr,eps,aeps)

            if (verbosity>=10 .and. depth<999) then
                print *,"recursive",depth
                print *,size(ml%weights)
                print *,tr(:,1)
                print *,tr(:,2)
                print *,score
                print *,"    ",best_score
                print *,"    ",best_tr(:,1)
                print *,"    ",best_tr(:,2)
                print *,deltav(tr)
            end if

            if (score(2)<minq) then
                print *,"prune",score(2),minq
                return
            end if

            if (all(deltav(tr)<mindelta)) then
                if (score(2)>best_score(2)) then
                    best_score = score
                    best_tr = tr
                end if
                return
            end if

            call decimate(eml,ml)

            call split(subtr,tr)

            do i=1,size(subtr,1)
                call recursive(eml,subtr(i,:,:),depth+1)
            end do
        end subroutine recursive
    end subroutine search_pp
end module rast
