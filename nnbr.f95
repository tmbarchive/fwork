!!! is there any way to do this per class?

module nnbr

    use dataio
    use quicksort

    implicit none

    integer :: verbose = 0
    real, parameter :: inf = 1e30
    real, parameter :: mindist = 1e-5

contains

    ! Given a list of values in x(:) , this finds the k-best values of that list
    ! (k is determined by the size of the output array).  This is the equivalent
    ! of the first k elements of argsort.

    subroutine kbest(neighbors,x)
        integer neighbors(:)
        real x(:)
        real vals(size(neighbors,1)),val
        integer k,n,i,j
        n = size(x,1)
        k = size(neighbors,1)
        if (k>n) stop "range error"
        neighbors = -1
        vals = inf
        do i=1,n
            val = x(i)
            j=k
            if (val<mindist) cycle
            if (val>=vals(k)) cycle
            do while (j>=1)
                if (j==1) then
                    vals(j) = val
                    neighbors(j) = i
                    exit
                else if (vals(j-1)<val) then
                    vals(j) = val
                    neighbors(j) = i
                    exit
                else
                    vals(j) = vals(j-1)
                    neighbors(j) = neighbors(j-1)
                    j = j-1
                end if
            end do
        end do
    end subroutine kbest

    ! compute the indexes of the k nearest neighbors from
    ! list of prototypes

    subroutine k_nearest_neighbors(is,x,protos)
        integer :: is(:)
        real :: x(:),protos(:,:)
        real :: dists(size(protos,1))
        integer i
        !$omp parallel do shared(dists,x,protos)
        do i=1,size(dists)
            dists(i) = sqrt(sum((x-protos(i,:))**2))
        end do
        !$omp end parallel do
        call kbest(is,dists)
        if (any(is>size(protos,1))) stop "oops"
    end subroutine k_nearest_neighbors

    ! compute the indexes of the k nearest neighbors, skipping
    ! prototypes in the skip list

    subroutine k_nearest_neighbors_skip(is,x,protos,skip)
        integer :: is(:)
        real :: x(:),protos(:,:)
        real :: dists(size(protos,1))
        logical :: skip(:)
        integer i
        !$omp parallel do shared(dists,x,protos)
        do i=1,size(dists)
            if(.not. skip(i)) then
                dists(i) = sqrt(sum((x-protos(i,:))**2))
            else
                dists(i) = 1e30
            end if
        end do
        !$omp end parallel do
        call kbest(is,dists)
        if (any(is>size(protos,1))) stop "oops"
    end subroutine k_nearest_neighbors_skip

    ! If nbrcls is a list of classes associated with some point x,
    ! compute the knn classifications for all k

    subroutine compute_knn_classes(classifications,nbrcls)
        integer c,nc,nbrcls(:),classifications(:),i,pred,actual,nclasses
        integer hist(maxval(nbrcls))
        classifications = 0
        hist = 0
        nclasses = maxval(nbrcls)
        do i=1,size(nbrcls)
            if (nbrcls(i)<1) cycle
            hist(nbrcls(i)) = hist(nbrcls(i))+1
            classifications(i) = maxloc(hist,1)
        end do
    end subroutine compute_knn_classes

    subroutine compute_knn_errors(errvec,nbrcls,actual)
        integer errvec(:),nbrcls(:),actual
        call compute_knn_classes(errvec,nbrcls)
        where (errvec==actual)
            errvec = 0
        elsewhere
            errvec = 1
        end where
    end subroutine compute_knn_errors

    integer function knn_best_k(protos,classes)
        integer, parameter :: kmax = 100
        real protos(:,:)
        integer classes(:)
        integer errs(kmax),errvec(kmax),neighbors(kmax)
        integer i
        errs = 0
        do i=1,size(protos,1)
            call k_nearest_neighbors(neighbors,protos(i,:),protos)
            neighbors = classes(neighbors)
            call compute_knn_errors(errvec,neighbors,classes(i))
            errs = errs + errvec
        end do
        knn_best_k = minloc(errs,1)
    end function knn_best_k


    ! Fill the picks(:) array sequentially with elements from list(:)
    ! chosen according to proability p

    ! TODO: rewrite this as COMPRESS(... random<=p ...)

    subroutine pick(picks,list,p)
        integer picks(:),list(:),i,j
        real p,q
        picks = -1
        i = 1
        j = 1
        do while (i<=size(picks) .and. j<=size(list))
            call random_number(q)
            if (q<=p) then
                picks(i) = list(j)
                i = i+1
            end if
            j = j+1
        end do
    end subroutine pick

    ! Evaluate k-nearest-neighbor errors for a single target vector x for a range of
    ! training set sizes and a range of k.  The results are returned in the accumulator
    ! array.  It contains the number of errors at different training set sizes n and
    ! different k in accumulator(n,k,1), and the corresponding counts in (:,:,2).
    ! The conditions array contains a range of sizes n to be tested in (:,1) and the
    ! corresponding number of repetitions in (:,2).

    subroutine evaluate_knn_errors(accumulator,conditions,x,actual,protos,classes)
        integer accumulator(:,:,:) ! (:,:,1) is errs, (:,:,2) is counts
        integer conditions(size(accumulator,1),2) ! (:,1) is n, (:,2) is repeat
        real x(:),protos(:,:),dists(size(protos,1)),p
        integer best(size(protos,1)),picks(size(accumulator,2)),k
        integer nbrcls(size(accumulator,2)),actual,classes(:),i,repeat
        integer nconditions,nprotos,n,r,errvec(size(accumulator,2)),nclasses
        k = size(accumulator,2)
        nclasses = maxval(classes)
        nconditions = size(conditions,1)
        nprotos = size(protos,1)

        !This loop dominates the entire computation; parallelizing it
        !results in linear speedups in the number of CPUs.  No other
        !part of this program matters in comparison.

        !$omp parallel shared(dists)
        !$omp do schedule(runtime)
        do i=1,size(dists)
            dists(i) = sqrt(sum((x-protos(i,:))**2))
        end do
        !$omp end do
        !$omp end parallel

        where(dists<mindist) dists = inf

        call quicksort1(best,dists,size(dists))

        if (size(conditions,2)/=2) stop "bad conditions"
        do i=1,nconditions
            n = conditions(i,1)
            repeat = conditions(i,2)
            p = n/real(nprotos)
            if (verbose>=100) print *,"condition",n,repeat,p
            if (n<nclasses) stop "bad n"
            if (repeat<1 .or. repeat>1000000) stop "bad repeat"
            do r=1,repeat
                call pick(picks,best,p)
                if (verbose>=100) print *,"trial",r,picks
                nbrcls = -1
                where (picks>0) nbrcls = classes(picks)
                call compute_knn_errors(errvec,nbrcls,actual)
                where (picks>0) accumulator(i,:,1) = accumulator(i,:,1) + errvec
                where (picks>0) accumulator(i,:,2) = accumulator(i,:,2) + 1
            end do
        end do
    end subroutine evaluate_knn_errors

    subroutine evaluate_knn_and_n(accumulator,conditions,protos,classes,maxtrials)
        real :: protos(:,:)
        integer :: classes(:)
        integer :: conditions(:,:)
        integer :: accumulator(:,:,:)
        integer, optional :: maxtrials
        integer nconditions,i,j,k
        real rate(size(accumulator,1),size(accumulator,2))
        integer ntrials

        nconditions = size(accumulator,1)
        if (nconditions/=size(conditions,1)) stop "accumulator and conditions inconsistent"
        k = size(accumulator,2)
        accumulator = 0
        ntrials = size(protos,1)
        if (present(maxtrials)) ntrials = min(ntrials,maxtrials)
        do i=1,ntrials
            if (verbose>=0 .and. modulo(i-1,100)==0) then
                write (*,fmt="(i0,1x)",advance="no") i
            end if
            call evaluate_knn_errors(accumulator,conditions,protos(i,:),classes(i),protos,classes)
            if (verbose>=10 .and. modulo(i,1000)==0) then
                rate = -1
                where (accumulator(:,:,2)>0) rate = accumulator(:,:,1) / real(accumulator(:,:,2))
                do j=1,size(rate,1)
                    print *,rate(j,1:min(k,8))
                end do
            end if
        end do
        if (verbose>=0) then
            write (*,fmt="(a)") "done"
        end if
    end subroutine evaluate_knn_and_n

    subroutine print_knn_and_n(accumulator,conditions)
        integer :: conditions(:,:)
        integer :: accumulator(:,:,:)
        real rate(size(accumulator,1),size(accumulator,2))
        integer j
        rate = -1
        where (accumulator(:,:,2)>0) rate = accumulator(:,:,1) / real(accumulator(:,:,2))
        do j=1,size(rate,1)
            write(*,fmt="(i7,a,i2,a)",advance="no") conditions(j,1)," [",conditions(j,2),"] "
            print '(99f8.4)',rate(j,:)
        end do
    end subroutine print_knn_and_n

end module nnbr
