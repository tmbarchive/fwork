program uni_mnist
    use dataio
    use uniclasses
    implicit none

    type(uniclass) uc
    real, allocatable :: inputs(:,:),classes(:)
    real, allocatable :: posterior(:)
    integer nfeat,nsamples,nclasses,i

    call read_mnist(inputs,"mnist/train-images-idx3-ubyte")
    call read_mnist(classes,"mnist/train-labels-idx1-ubyte")

    inputs = inputs / maxval(inputs)
    nfeat = size(inputs,2)
    nsamples = size(classes)
    classes = classes + 1
    nclasses = maxval(classes)
    allocate(posterior(nclasses))

    call uniclass_init(uc,nfeat,nclasses)

    do i=1,nsamples
        if(modulo(i,1000)==0) then
            call uniclass_posterior(uc,posterior,inputs(i,:))
            print *,i,floor(classes(i)-1),maxloc(posterior)-1
        end if
        call uniclass_add(uc,inputs(i,:),floor(0.5+classes(i)))
    end do
    call uniclass_retrain(uc)
    call uniclass_print(uc,1)
end program uni_mnist
