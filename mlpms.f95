module mlps

    use dataio
    use quicksort
    use mlps

    implicit none

    type mlpmulti
        integer :: nnets = 0
        type(mlp), allocatable :: nets(:)
    end type mlpmulti

contains

    subroutine mlpm_write(nets,u)
        type(mlpmulti) nets
        integer u
        do i=1,nets%nnets
            call mlp_write(nets%nets(i),u)
        end do
    end subroutine mlpm_write

    subroutine mlpm_read(nets,u)
        type(mlpmulti) nets
        integer u
        do i=1,nets%nnets
            call mlp_read(nets%nets(i),u)
        end do
    end subroutine mlpm_read

    subroutine mlpm_forward(nets,z,x)
        type(mlp) net
        real, intent(out) :: z(:)
        real, intent(in) :: x(:)
        real results(nets%nnets,size(z))
        !$omp parallel do shared(results,z,x)
        do i=1,nets%nnets
            call mlp_forward(nets%nets(i),results(i,:),x)
        end do
        !$omp end parallel do
        z = sum(results,1) / size(results,1)
    end subroutine mlpm_forward

    subroutine mlpm_train1(nets,targets,inputs,frac)
        type(autoparams) ps
        type(mlpmulti) nets
        real, intent(in) :: targets(:,:)
        real, intent(in) :: inputs(:,:)
        real :: frac
        if (.not. allocated(nets%nets)) then
            nets%nnets = 0
            allocate(nets%nets(100))
        end if
        nets%nnets = nets%nnets + 1
        call mlp_autotrain(nets%nets(nets%nnets),targets,inputs,frac)
    end subroutine mlpm_train1

    subroutine mlpm_train_cv(nets,targets,inputs,frac,n)
        type(autoparams) ps
        type(mlpmulti) nets
        real, intent(in) :: targets(:,:)
        real, intent(in) :: inputs(:,:)
        real :: frac
        integer n
        if (allocated(nets%nets)) deallocate(nets%nets)
        do i=1,n
            call mlp_autotrain(nets%nets(nets%nnets),targets,inputs,frac)
        end do
    end subroutine mlpm_train_cv

end module mlpms
