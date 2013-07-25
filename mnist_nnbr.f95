!!!
!!! TODO:
!!! -- rescaling
!!! -- normalization
!!! -- feature extractors
!!! -- edited nearest neighbor
!!! -- extrapolation to determine asymptotic error;
!!!    error estimate via resampling
!!! -- optimal stopping...

program mnist_nnbr
    use dataio
    use nnbr
    use mlps
    use fextract
    implicit none
    type(mlp) net
    real, allocatable :: images(:,:,:),iclasses(:),classes(:)
    real, allocatable :: vectors(:,:)
    integer nsamples,nclasses,nfeat
    integer epoch,i,which,ndim,n,j
    integer, parameter :: kmax = 21
    integer, parameter :: nconditions = 6
    integer, parameter :: maxtrials = 5000
    integer conditions(nconditions,2)
    integer accumulator(nconditions,kmax,2)
    integer, allocatable :: perm(:)
    real v(10000)


    ! load the mnist data
    call read_mnist(images,"mnist/train-images-idx3-ubyte")
    call read_mnist(iclasses,"mnist/train-labels-idx1-ubyte")
    nsamples = size(iclasses)
    iclasses = iclasses + 1
    nclasses = maxval(iclasses)
    if (size(images,1)/=size(iclasses)) stop "oops nsamples"
    images = images / maxval(images)

    ! generate a permutation
    allocate(perm(nsamples))
    allocate(classes(nsamples))
    call permutation1(perm,nsamples)


    conditions = -999999
    conditions(:,1) = [10000,20000,30000,40000,50000,60000]
    conditions(:,2) = [10,10,10,10,10,10]

    do which=1,8
        print *,"extract",which
        call extract(v,ndim,images(1,:,:),which)
        allocate(vectors(nsamples,ndim))
        do i=1,nsamples
            j = perm(i)
            call extract(vectors(i,:),n,images(j,:,:),which)
            classes(i) = iclasses(j)
            if (n/=ndim) stop "oops ndim"
        end do
        print *,"results",which
        accumulator = -1
        call evaluate_knn_and_n(accumulator,conditions,vectors,floor(classes),&
            &maxtrials=maxtrials)
        call print_knn_and_n(accumulator,conditions)
        deallocate(vectors)
    end do
end program mnist_nnbr
