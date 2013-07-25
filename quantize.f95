module quantize

    use utils
    use quicksort

    implicit none

    ! A quantizer used in a way similar to a nearest
    ! neighbor classifier; the main difference is that
    ! we store posteriors rather than classes associated
    ! with each prototype.  ntraining*posteriors is the
    ! counts that went into each prototype.

    type quantizer
        integer :: ntraining
        real, allocatable :: protos(:,:)
        real, allocatable :: posteriors(:,:)
        real :: rate1
        real :: rate2
    end type quantizer

    type quantparams
        integer :: nquants = 4
        integer :: stagesize = 20000
        integer :: nstages = 4
    end type quantparams

    integer :: verbose = 0

contains

    ! standard k-means algorithm; too slow for large
    ! datasets

    subroutine kmeans(means,data,maxiter,k,n,d)
        integer, intent(in) :: maxiter,k,n,d
        real, intent(out) :: means(k,d)
        real, intent(in) ::  data(n,d)
        integer assignment(n),last(n),counts(k)
        integer i,c,iter
        print *,"kmeans",k,n,d
        do i=1,k
            means(i,:) = data(i*(n/k),:)
        end do
        last(:) = 0
        do iter=1,maxiter
            print *,"iter",iter
            do i=1,n
                assignment(i) = argmindist(data(i,:),means)
            end do
            counts(:) = 0
            means(:,:) = 0
            do i=1,n
                c = assignment(i)
                counts(c) = counts(c) + 1
                means(c,:) = means(c,:) + data(i,:)
            end do
            if (all(assignment(:)==last(:))) exit
            do i=1,k
                if (counts(i)/=0) then
                    means(i,:) = means(i,:) / real(counts(i))
                else
                    means(i,:) = data(i,:) ! bad initializer
                end if
            end do
            last(:) = assignment(:)
        end do
        print *,"done"
    end subroutine kmeans

    ! incremental kmeans algorithm
    ! TODO: speed up with bounds

    subroutine incremental_kmeans(means,data,maxiter,rate,k,n,d)
        integer, intent(in) :: maxiter,k,n,d
        real, intent(out) :: means(k,d)
        real, intent(in) ::  data(n,d), rate
        integer c,i2,j,iter
        real d1,d2,dist,l
        do c=1,k
            j = randint(n)
            means(c,:) = data(j,:)
        end do
        dist = 0.0
        l = 1.0/1024.0
        do iter=1,maxiter
            if (verbose>=10 .and. modulo(iter,5000)==0) print *,"ikm",iter,dist
            j = randint(n)
            call argmindist2(c,d1,i2,d2,data(j,:),means,k,d)
            means(c,:) = means(c,:) + rate * (data(j,:) - means(c,:))
            dist = (1.0-l) * dist + l * d1
        end do
    end subroutine incremental_kmeans

    ! neural gas algorithm; basically the same as kmeans, but
    ! has a little bit of "regularization" in it (may lead
    ! to better balanced sets, fewer outliers)

    ! TODO: speed up with bounds

    subroutine neural_gas_init(means,data,k)
        real, allocatable, intent(out) :: means(:,:)
        real, intent(in) ::  data(:,:)
        integer :: maxiter,n,k,i,i1,d
        n = size(data,1)
        d = size(data,2)
        allocate(means(k,d))
        do i=1,k
            i1 = randint(n)
            means(i,:) = data(i1,:)
        end do
    end subroutine neural_gas_init

    subroutine neural_gas(means,data,maxiter,rate1,rate2)
        integer, intent(in) :: maxiter
        real, intent(out) :: means(:,:)
        real, intent(in) ::  data(:,:), rate1, rate2
        integer i1,i2,iter,i,j,k,n,d
        real d1,d2,dist,l
        k = size(means,1)
        n = size(data,2)
        d = size(means,2)
        dist = 0.0
        l = 1.0/1024.0
        do iter=1,maxiter
            if (verbose>=10 .and. modulo(iter,10000)==0) print *,"ngq",iter,dist
            j = randint(n)
            call argmindist2(i1,d1,i2,d2,data(j,:),means,k,d)
            means(i1,:) = means(i1,:) + rate1 * (data(j,:) - means(i1,:))
            means(i2,:) = means(i2,:) + rate2 * (data(j,:) - means(i2,:))
            dist = (1.0-l) * dist + l * d1
        end do
    end subroutine neural_gas

    function quant_eval(quant,data) result(q)
        type(quantizer) quant
        integer :: k,n,i
        real :: data(:,:)
        real :: dists(size(data,1))
        real q
        k = size(quant%protos,1)
        n = size(data,1)
        do i=1,n
            dists(i) = mindist(data(i,:),quant%protos)
        end do
        q = sum(dists)/size(dists)
    end function quant_eval

    function quant_eval_impurity(quant,data,classes) result(q)
        type(quantizer) quant
        real :: data(:,:)
        integer classes(:)
        integer :: k,n,i
        real :: dists(size(quant%protos,1))
        real q,d
        stop 'unimplemented'
    end function quant_eval_impurity

    subroutine quant_train(qs,quant,data,k)
        type(quantparams) qs
        type(quantizer) quant,quants(qs%nquants)
        real :: rates(qs%nquants,2),errs(qs%nquants)
        real :: data(:,:),best
        integer i,j,k,stage,order(qs%nquants)

        do j=1,qs%nquants
            call neural_gas_init(quants(j)%protos,data,k)
            quants(j)%rate1 = exp(randnormal()+log(0.1))
            quants(j)%rate2 = exp(randnormal()+log(0.02))
            quants(j)%rate2 = min(quants(j)%rate1,quants(j)%rate2)
        end do

        if(allocated(quant%protos)) deallocate(quant%protos)
        if(allocated(quant%posteriors)) deallocate(quant%posteriors)

        best = huge(1.0)
        do stage=1,qs%nstages
            errs = -1
            do j=1,qs%nquants
                call neural_gas(quants(j)%protos,data,qs%stagesize,quants(j)%rate1,quants(j)%rate2)
                errs(j) = quant_eval(quants(j),data)
                print *,errs
            end do
            call quicksort1(order,errs,size(errs))
            print *,"stage",stage
            print *,"errs",(errs(order(j)),j=1,qs%nquants)
            print *,"etas",(quants(order(j))%rate1,j=1,qs%nquants)
            print *,"etas",(quants(order(j))%rate2,j=1,qs%nquants)
            if (errs(order(1))<best) then
                best = errs(order(1))
                quant = quants(order(1)) ! save the best net so far
            end if
            do i=1,qs%nquants/2
                j = order(i+qs%nquants/2)
                quants(j) = quants(order(i))
                quants(j)%rate1 = exp(randnormal()+log(quants(j)%rate1))
                quants(j)%rate2 = exp(randnormal()+log(quants(j)%rate2))
                quants(j)%rate2 = min(quants(j)%rate1,quants(j)%rate2)
            end do
        end do

    end subroutine quant_train

    subroutine quant_train_impurity(qs,quant,data,k,classes)
        type(quantparams) qs
        type(quantizer) quant,quants(qs%nquants)
        integer :: classes(:)
        real :: rates(qs%nquants,2),errs(qs%nquants)
        real :: data(:,:),best
        integer i,j,k,stage,order(qs%nquants)

        do j=1,qs%nquants
            call neural_gas_init(quants(j)%protos,data,k)
            quants(j)%rate1 = exp(randnormal()+log(0.1))
            quants(j)%rate2 = exp(randnormal()+log(0.02))
            quants(j)%rate2 = min(quants(j)%rate1,quants(j)%rate2)
        end do

        if(allocated(quant%protos)) deallocate(quant%protos)
        if(allocated(quant%posteriors)) deallocate(quant%posteriors)

        best = huge(1.0)
        do stage=1,qs%nstages
            errs = -1
            do j=1,qs%nquants
                call neural_gas(quants(j)%protos,data,qs%stagesize,quants(j)%rate1,quants(j)%rate2)
                errs(j) = quant_eval_impurity(quants(j),data,classes)
                print *,errs
            end do
            call quicksort1(order,errs,size(errs))
            print *,"stage",stage
            print *,"errs",(errs(order(j)),j=1,qs%nquants)
            print *,"etas",(quants(order(j))%rate1,j=1,qs%nquants)
            print *,"etas",(quants(order(j))%rate2,j=1,qs%nquants)
            if (errs(order(1))<best) then
                best = errs(order(1))
                quant = quants(order(1)) ! save the best net so far
            end if
            do i=1,qs%nquants/2
                j = order(i+qs%nquants/2)
                quants(j) = quants(order(i))
                quants(j)%rate1 = exp(randnormal()+log(quants(j)%rate1))
                quants(j)%rate2 = exp(randnormal()+log(quants(j)%rate2))
                quants(j)%rate2 = min(quants(j)%rate1,quants(j)%rate2)
            end do
        end do

    end subroutine quant_train_impurity

end module quantize
