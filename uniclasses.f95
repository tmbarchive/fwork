!!!
!!! uniform classifiers
!!! just add samples, and it picks the right representation internally
!!!
!!! TODO:
!!! add retrain method
!!! add feature/preprocessing selection
!!! compute knn error, use k for classification, error for bp limit

module uniclasses
    use dataio
    use nnbr
    use knncs
    use mlps
    use quantize
    use utils

    implicit none

    type uniclass
        character(1000) id
        integer nfeat, nclasses
        type(knnc), pointer :: knn
        type(mlp), pointer :: net => null()
        type(quantizer), pointer :: quant => null()
        type(uniclass), pointer :: sub(:) => null()
    end type uniclass

    integer :: mlp_min = 2000
    integer :: mlp_threshold = 100000
    integer :: split_threshold = 50000
    integer :: nsplit = 10

    logical :: initialized = .false.

    type(autoparams), save :: retrain_ps
    type(autoparams), save :: online_ps
    type(quantparams), save :: split_qs
contains

    subroutine uniclass_globals
        retrain_ps%stagesize = 100000
        retrain_ps%nstages = 10
        retrain_ps%nnets = 4
        retrain_ps%hidden_lo = 5
        retrain_ps%hidden_hi = 50

        online_ps%stagesize = 50000
        online_ps%nstages = 10
        online_ps%nnets = 4
        online_ps%hidden_lo = 10
        online_ps%hidden_hi = 50

        split_qs%nquants = 4
        split_qs%stagesize = 10000
        split_qs%nstages = 3
    end subroutine uniclass_globals

    subroutine uniclass_check(uc)
        type(uniclass) uc
        integer n
        if (uc%nclasses<2 .or. uc%nclasses>1000) stop "oops nclasses"
        if (uc%nfeat<2 .or. uc%nfeat>1000000) stop "oops nfeat"

        if (associated(uc%knn)) then
            if (uc%knn%nprotos>0) then
                n = uc%knn%nprotos
                if (uc%nfeat/=size(uc%knn%protos,2)) stop "oops nfeat"
                if (minval(uc%knn%classes(:n))<1) stop "oops knn classes lo"
                if (maxval(uc%knn%classes(:n))>uc%nclasses) stop "oops knn classes hi"
                if (size(uc%knn%protos(:n,:),1)<1) stop "oops"
                if (size(uc%knn%protos(:n,:),1)>1000000) stop "oops"
                if (maxval(abs(uc%knn%protos(:n,:)))>100) stop "oops"
            end if
        else
            stop "oops knn not associated"
        end if

        if (associated(uc%net)) then
            if (uc%nfeat/=uc%net%ninput) stop "oops mlp"
            if (uc%nclasses/=uc%net%noutput) stop "oops mlp"
        end if

        if (associated(uc%quant)) then
            if (uc%nfeat/=size(uc%quant%protos,2)) stop "oops nfeat"
            if (size(uc%quant%protos,1)<1) stop "oops"
            if (size(uc%quant%protos,1)>1000) stop "oops"
            if (maxval(abs(uc%quant%protos))>100) stop "oops"
        end if

        if (associated(uc%sub)) then
            if (size(uc%sub)/=size(uc%quant%protos,1)) stop "oops nsub"
        end if
    end subroutine uniclass_check

    subroutine uniclass_init(uc,nfeat,nclasses)
        type(uniclass) uc
        integer nfeat,nclasses
        uc%id = "top"
        uc%nfeat = nfeat
        uc%nclasses = nclasses
        allocate(uc%knn)
    end subroutine uniclass_init

    recursive subroutine uniclass_retrain(uc)
        type(uniclass) uc
        integer i,n
        real, allocatable :: cm(:,:)
        integer, allocatable :: h(:)

        n = uc%knn%nprotos
        print *,"==============="
        print *,"retraining ",trim(uc%id),n
        print *,"==============="
        if (associated(uc%sub)) then
            do i=1,size(uc%sub)
                call uniclass_retrain(uc%sub(i))
            end do
        else
            if (n>=mlp_min) then
                call knn_seal(uc%knn)
                call print_hist(uc%knn%classes)
                call hist(h,uc%knn%classes)

                if(count(h>0)>1) then
                    if (associated(uc%net)) deallocate(uc%net)
                    allocate(uc%net)
                    call mlp_autotrain_cv(retrain_ps,uc%net,uc%knn%classes,uc%knn%protos)
                    call mlp_check(uc%net)
                    call mlp_confusion(cm,uc%net,uc%knn%classes,uc%knn%protos)
                    call print_confusion(cm)
                end if
            end if
        end if
    end subroutine uniclass_retrain

    recursive subroutine uniclass_print(uc,depth)
        type(uniclass) uc
        integer depth,i
        integer, allocatable :: h(:)
        call indent(depth)
        call hist(h,uc%knn%classes)
        print *,uc%knn%nprotos,h
        if (associated(uc%sub)) then
            do i=1,size(uc%sub)
                call uniclass_print(uc%sub(i),depth+1)
            end do
        end if
    end subroutine uniclass_print

    subroutine uniclass_split(uc)
        type(uniclass) uc
        integer c,n,bucket,i
        character(1000) buffer

        call uniclass_check(uc)
        print *,"==============="
        print *,"uniclass splitting: ",trim(uc%id)
        print *,"==============="
        n = uc%knn%nprotos
        call knn_seal(uc%knn)
        call print_hist(uc%knn%classes)
        allocate(uc%quant)
        call quant_train(split_qs,uc%quant,uc%knn%protos,nsplit)
        allocate(uc%sub(nsplit))
        do i=1,nsplit
            call uniclass_init(uc%sub(i),uc%nfeat,uc%nclasses)
            write (uc%sub(i)%id,fmt="(a,a,i2)") trim(uc%id),"/",i
        end do
        do i=1,n
            bucket = argmindist(uc%knn%protos(i,:),uc%quant%protos)
            call uniclass_add(uc%sub(bucket),uc%knn%protos(i,:),uc%knn%classes(i))
        end do
        do i=1,nsplit
            call uniclass_check(uc)
            n = uc%sub(i)%knn%nprotos
            if (n<1) cycle
            print *,"*** bucket",i,size(uc%sub(i)%knn%classes(:n)),&
                &maxval(uc%sub(i)%knn%classes(:n))
            ! call print_hist(uc%sub(i)%knn%classes(:n))
        end do
    end subroutine uniclass_split

    recursive subroutine uniclass_add(uc,v,c)
        character(1000) buffer
        type(uniclass) uc
        real v(:)
        integer c,n,bucket,i

        if (.not.initialized) then
            initialized = .true.
            call uniclass_globals
        end if
        if (c<1 .or. c>9999) stop "bad class"
        if (maxval(abs(v))>100) stop "bad feature vector"

        if (associated(uc%sub)) then
            bucket = argmindist(v,uc%quant%protos)
            call uniclass_add(uc%sub(bucket),v,c)
        else

            if (associated(uc%quant)) then
                bucket = argmindist(v,uc%quant%protos)
                call uniclass_add(uc%sub(bucket),v,c)
            else
                call knn_add(uc%knn,v,c)
                n = uc%knn%nprotos
                if (n>=split_threshold) then
                    call uniclass_split(uc)
                end if
            end if
        end if
    end subroutine uniclass_add

    recursive subroutine uniclass_posterior(uc,posterior,v)
        type(uniclass) uc
        real posterior(:)
        real v(:)
        integer bucket,n

        if (.not.associated(uc%net) .and. uc%knn%nprotos>=mlp_threshold) then
            ! just do a quick and dirty job on the training
            print *,"==============="
            print *,"uniclass training mlp: ",trim(uc%id)
            print *,"==============="
            allocate(uc%net)
            n = uc%knn%nprotos
            call mlp_autotrain_cv(online_ps,uc%net,uc%knn%classes(:n),uc%knn%protos(:n,:))
            call mlp_check(uc%net)
        end if

        if (associated(uc%sub)) then
            bucket = argmindist(v,uc%quant%protos)
            call uniclass_posterior(uc%sub(bucket),posterior,v)
        else if (associated(uc%net)) then
            call mlp_forward(uc%net,posterior,v)
        else
            call knn_posterior(uc%knn,posterior,v)
        end if
    end subroutine uniclass_posterior

    function uni_error(uc,classes,inputs) result(errs)
        type(uniclass) uc
        integer, intent(in) :: classes(:)
        real, intent(in) :: inputs(:,:)
        real errs

        errs = 0
    end function uni_error

end module uniclasses
