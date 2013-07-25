program mlp_mnist
    use dataio
    use mlps
    type(mlp) net
    real, allocatable :: inputs(:,:),classes(:)
    integer nsamples,nclasses,nfeat
    integer epoch,i

    call read_mnist(inputs,"mnist/train-images-idx3-ubyte")
    call read_mnist(classes,"mnist/train-labels-idx1-ubyte")

    nsamples = size(classes)
    classes = classes + 1
    nclasses = maxval(classes)
    if (size(images,1)/=size(classes)) stop "oops nsamples"
    images = images / maxval(images)

    do which=5,8
        print *,"extract",which
        call extract(v,ndim,images(1,:,:),which)
        allocate(vectors(nsamples,ndim))
        do i=1,nsamples
            call extract(vectors(i,:),n,images(i,:,:),which)
            if (n/=ndim) stop "oops ndim"
        end do

!!$        call mlp_autotrain here
!!$        call mlp_init(net,ndim,20,nclasses)
!!$        do epoch=1,100
!!$            print *,"***",epoch
!!$            call mlp_clear_info(net)
!!$            call mlp_train(net,reshape(classes,[nsamples,1]),inputs,0.3)
!!$            print *,epoch,mlp_error(net,reshape(classes(1:1000),[1000,1]),inputs(1:1000,:))
!!$            call mlp_print_info(net)
!!$        end do
        deallocate(vectors)
    end do
end program mlp_mnist
