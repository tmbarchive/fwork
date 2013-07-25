module knncs
    use reallocs
    use nnbr

    implicit none

    ! this data type serves both as a classifier
    ! and as a repository for training data for other
    ! classifiers

    type knnc
        integer :: nprotos = 0
        integer :: k = 0
        real, allocatable :: protos(:,:)
        integer, allocatable :: classes(:)
    end type knnc

    integer, parameter :: initial_allocation = 10000

contains

    ! Given a list of elements, find the most frequent element.  This is
    ! for direct computation of nearest neighbor classifiers.

    integer function most_frequent(classes)
        integer classes(:),hist(1000),nc,c,i
        nc = maxval(classes)
        hist(1:nc) = 0
        do i=1,size(classes)
            c = classes(i)
            hist(c) = hist(c)+1
        end do
        most_frequent = maxloc(hist,1)
    end function most_frequent

    ! Compute the distance of x from each of the rows of the protos array,
    ! then compute the indexes of the k closest rows  and return it.

    ! estimate the best k for knn classification

    subroutine knn_estimate_k(knn,nsamples)
        type(knnc) knn
        integer, optional :: nsamples
        integer n
        n = knn%nprotos
        if (present(nsamples)) n = min(knn%nprotos,nsamples)
        knn%k = knn_best_k(knn%protos(1:n,:),knn%classes(1:n))
    end subroutine knn_estimate_k

    ! "train" the nearest neighbor classifier (just stores the vector)

    subroutine knn_add(knn,v,c)
        type(knnc) knn
        real v(:)
        integer c
        if (.not. allocated(knn%protos)) then
            knn%nprotos = 0
            knn%k = 1
            allocate(knn%protos(initial_allocation,size(v)))
            knn%protos = huge(1.0)
            allocate(knn%classes(initial_allocation))
            knn%classes = huge(1)
        else if (knn%nprotos==size(knn%protos,1)) then
            print *,"[ realloc",knn%nprotos,"]"
            call realloc(knn%protos,knn%nprotos*2,size(knn%protos,2))
            call realloc(knn%classes,knn%nprotos*2)
        end if
        knn%nprotos = knn%nprotos+1
        knn%protos(knn%nprotos,:) = v
        knn%classes(knn%nprotos) = c
    end subroutine knn_add

    ! reallocate everything inside the knn structure to have the right size

    subroutine knn_seal(knn)
        type(knnc) knn
        call realloc(knn%protos,knn%nprotos,size(knn%protos,2))
        call realloc(knn%classes,knn%nprotos)
    end subroutine knn_seal

    ! classify with the nearest neighbor classifier

    integer function knn_classify(knn,v) result(c)
        type(knnc) knn
        integer neighbors(knn%k)
        real v(:)
        call k_nearest_neighbors(neighbors,v,knn%protos)
        neighbors = knn%classes(neighbors)
        c = most_frequent(neighbors)
    end function knn_classify

    ! compute an estimate of the posterior with the nearest
    ! neighbor classifier

    subroutine knn_posterior(knn,posterior,v)
        type(knnc) knn
        real v(:),posterior(:)
        integer neighbors(knn%k),cls(knn%k)
        integer i,k,nprotos
        nprotos = knn%nprotos
        k = knn%k
        posterior = 0
        call k_nearest_neighbors(neighbors,v,knn%protos(:nprotos,:))
        if (any(neighbors>nprotos)) stop "oops"
        do i=1,k
            cls(i) = knn%classes(neighbors(i))
        end do
        do i=1,k
            posterior(cls(i)) = posterior(cls(i))+1
        end do
        posterior = posterior / sum(posterior)
    end subroutine knn_posterior
end module knncs
