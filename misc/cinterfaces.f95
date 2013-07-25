!!!
!!! C-callable interfaces for the classifiers defined in this package.
!!!

module fclassifier_globals
    use mlps
    use uniclasses
    implicit none
    type(mlp) nets(0:10)
    type(uniclass) ucs(0:10)
end module fclassifier_globals

module fclassifier
    use mlps
    use fclassifier_globals
    use iso_c_binding
    implicit none
    character(100) error
contains

    !!
    !! uniform classifier interface
    !!

    subroutine f_uniclass_add(uc,vector,cls,dim)
        integer(C_INT) uc,cls,dim
        real(C_FLOAT) vector(dim)
        call uniclass_add(ucs(uc),vector,cls)
    end subroutine f_uniclass_add

    subroutine f_uniclass_posterior(uc,posterior,vector,nposterior,dim)
        integer(C_INT) uc,nposterior,dim
        real(C_FLOAT) posterior(nposterior),vector(dim)
        call uniclass_posterior(ucs(uc),posterior,vector)
    end subroutine f_uniclass_posterior

    !!
    !! interfaces to Fortran MLP code
    !!
    !! NB: these routines expect class labels
    !! in the range 1...C (not 0...C-1), arrays in Fortran order
    !!

    subroutine f_mlp_init(net,ninput,nhidden,noutput) bind(c)
        integer(C_INT), value :: net,ninput,nhidden,noutput
        call mlp_init(nets(net),ninput,nhidden,noutput)
    end subroutine f_mlp_init

    subroutine f_mlp_train(net,targets,inputs,eta,nrows,ninputs,noutputs) bind(c)
        integer(C_INT), value :: net,nrows,ninputs,noutputs
        real(C_FLOAT), value :: eta
        real(C_FLOAT) targets(nrows,noutputs),inputs(nrows,ninputs)
        if(ninputs/=nets(net)%ninput) stop "wrong input size"
        if(noutputs/=1) then
            if(noutputs/=nets(net)%noutput) stop "wrong output size"
        else
            if(minval(targets)<1) stop "zero or negative class label"
            if(maxval(targets)>nets(net)%noutput) stop "class label too big"
        end if
        call mlp_train(nets(net),targets,inputs,eta)
    end subroutine f_mlp_train

    subroutine f_mlp_forward(net,z,x,noutputs,ninputs) bind(c)
        integer(C_INT), value :: net,noutputs,ninputs
        real(C_FLOAT) :: z(noutputs),x(ninputs)
        if(ninputs/=nets(net)%ninput) stop "wrong input size"
        if(noutputs/=nets(net)%noutput) stop "wrong output size"
        call mlp_forward(nets(net),z,x)
    end subroutine f_mlp_forward

    subroutine f_mlp_error(net,errs,targets,inputs,nrows,ninputs,noutputs) bind(c)
        integer(C_INT), value :: net,ninputs,noutputs,nrows
        real(C_FLOAT) :: targets(nrows,noutputs),inputs(nrows,ninputs),errs(2)
        if(ninputs/=nets(net)%ninput) stop "wrong input size"
        if(noutputs/=1) then
            if(noutputs/=nets(net)%noutput) stop "wrong output size"
        else
            if(minval(targets)<1) stop "zero or negative class label"
            if(maxval(targets)>nets(net)%noutput) stop "class label too big"
        end if
        errs = mlp_error(nets(net),targets,inputs)
    end subroutine f_mlp_error
end module fclassifier

