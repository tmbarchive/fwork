module fextract
    use cedges
    use dataio
    implicit none
contains

    subroutine extract1(v,n,im)
        real im(:,:)
        integer n
        real v(:)
        n = size(im)
        v(1:n) = reshape(im,[size(im)]) / sum(abs(im))
    end subroutine extract1

    subroutine extract2(v,n,im)
        real im(:,:)
        integer n
        real v(:)
        n = size(im)
        v(1:n) = reshape(im,[size(im)]) / sqrt(sum(im**2))
    end subroutine extract2

    subroutine extract3(v,n,im)
        real im(:,:)
        integer n
        real v(:)
        n = size(im)
        v(1:n) = reshape(im,[size(im)]) / sum(abs(im)**3)**(1.0/3)
    end subroutine extract3

    subroutine extractinf(v,n,im)
        real,intent(in) :: im(:,:)
        integer n
        real v(:)
        n = size(im)
        v(1:n) = reshape(im,[size(im)]) / maxval(im)
    end subroutine extractinf

    subroutine smooth1(v,n,im)
        real,intent(in) :: im(:,:)
        real :: image(size(im,1),size(im,2))
        integer n
        real v(:)
        n = size(im)
        image = im
        call gauss2d(image,1.0)
        ! call display(image)
        v(1:n) = reshape(image,[size(image)]) / sqrt(sum(image**2))
    end subroutine smooth1

    subroutine smooth2(v,n,im)
        real,intent(in) :: im(:,:)
        real :: image(size(im,1),size(im,2))
        integer n
        real v(:)
        n = size(im)
        image = im
        call gauss2d(image,2.0)
        v(1:n) = reshape(image,[size(image)]) / sqrt(sum(image**2))
    end subroutine smooth2

    subroutine smooth3(v,n,im)
        real,intent(in) :: im(:,:)
        real :: image(size(im,1),size(im,2))
        integer n
        real v(:)
        n = size(im)
        image = im
        call gauss2d(image,3.0)
        v(1:n) = reshape(image,[size(image)]) / sqrt(sum(image**2))
    end subroutine smooth3

    subroutine smooth4(v,n,im)
        real,intent(in) :: im(:,:)
        real :: image(size(im,1),size(im,2))
        integer n
        real v(:)
        n = size(im)
        image = im
        call gauss2d(image,4.0)
        v(1:n) = reshape(image,[size(image)]) / sqrt(sum(image**2))
    end subroutine smooth4

    subroutine extract(v,n,im,which)
        real, intent(in) :: im(:,:)
        real v(size(im))
        integer n
        integer which

        select case(which)
        case (1)
            call extract1(v,n,im)
        case (2)
            call extract2(v,n,im)
        case (3)
            call extract3(v,n,im)
        case (4)
            call extractinf(v,n,im)
        case (5)
            call smooth1(v,n,im)
        case (6)
            call smooth2(v,n,im)
        case (7)
            call smooth3(v,n,im)
        case (8)
            call smooth4(v,n,im)
        case default
            stop "bad extraction method"
        end select
    end subroutine extract

end module fextract

