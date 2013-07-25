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

    type pp_problem
        real eps,aeps
        type(point), allocatable :: model(:)
        type(point), allocatable :: image(:)
        integer(2), allocatable :: pairs(:,:)
        real, allocatable :: weights(:)
    end type pp_problem

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

    ! Evaluate point-point matches, given the selected list
    ! of pairs and the region in transformation space.
    
    function evaluate_pp(problem,list,tr) result(score)
        type(pp_problem) problem
        integer list(:)
        real tr(4,2),score(2)

        type(point) mpt,ipt
        integer last_mi,index,li,mi,ii
        real rfixed,rscaled,rtotal,tx,ty,tangle,tcos,tsin,tscale
        real mpx,mpy,mpa,d,d_scores(2),a_scores(2),total(2)

        ! compute the midpoints of the transformation
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

        ! iterate through all the model/image pairs and compute
        ! lower and upper bounds on their contributions to the score
        last_mi = -1
        total = 0
        do index=1,size(list)
            li = list(index)
            mi = problem%pairs(li,1)
            ii = problem%pairs(li,2)
            if (last_mi/=mi) then
                ! project the model point
                mpt = problem%model(mi)
                mpx = tcos * mpt%x - tsin * mpt%y + tx
                mpy = tsin * mpt%x + tcos * mpt%y + ty
                mpa = mpt%a + tangle
                last_mi = mi
                rtotal = rfixed + hypot(mpx,mpy) * rscaled
            end if
            ! todo: per-feature scores
            ipt = problem%image(ii)
            d = hypot(mpt%x-ipt%x,mpt%y-ipt%y)
            d_scores = d_score(d,rtotal,problem%eps) * problem%weights(li)
            if (d_scores(2)<1e-6) then
                list(index) = -1
            else
                total = total + d_scores
            end if
        end do
        score = total
    end function evaluate_pp

    subroutine search_pp(problem,list,tr,minq)
        real, parameter :: mindelta(4) = [1.0,1.0,0.01,0.002]
        type(pp_problem) problem
        integer list(:)
        real tr(4,2),minq
        real best_score(2),best_tr(4,2)

        if (problem%eps<0.1) stop "bad eps"
        if (problem%aeps<0.0001) stop "bad aeps"
        if (minval(problem%weights)<0.0) stop "bad weights"
        if (size(problem%weights)/=size(problem%pairs,1)) stop "mismatch"
        call search(list,tr,1)
    contains
        recursive subroutine search(list,tr,depth)
            integer list(:)
            integer sublist(count(list>0))
            real tr(4,2),score(2)
            real, allocatable :: subtr(:,:,:)
            integer depth,sub

            if (depth>100) stop "max depth exceeded"
            sublist = pack(list,list>0)
            score = evaluate_pp(problem,sublist,tr)
            if (depth<100) then
                print *,"===",depth,size(list)
                print *,"    ",score
                print *,"    ",tr(:,1)
                print *,"    ",tr(:,2)
                print *,"    best: ",best_score
            end if
            if (score(2)<minq) then
                print *,"pruned",score,minq
                print *,"    ",tr(:,1)
                print *,"    ",tr(:,2)
            end if
            if (score(2)<best_score(1)) then
                print *,"pruned2",score,minq
                print *,"    ",tr(:,1)
                print *,"    ",tr(:,2)
            end if
            if (all(deltav(tr)<mindelta)) then
                if (score(2)>best_score(2)) then
                    best_score = score
                    best_tr = tr
                    print *,"    ",tr(:,1)
                    print *,"    ",tr(:,2)
                end if
                return
            end if
            call split(subtr,tr)
            do sub=1,size(subtr,1)
                call search(sublist,subtr(sub,:,:),depth+1)
            end do
        end subroutine search
    end subroutine search_pp
end module rast
