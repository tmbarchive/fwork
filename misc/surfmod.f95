module surfmod
    type surfdesc
        real x,y,a,b,c,l
        real, allocatable :: descriptor(:)
    end type surfdesc
contains
    subroutine surf_read(d,fname)
        character(*) fname
        type(surfdesc), allocatable :: d(:)
        integer ndesc,npoints

        open(unit=10,file=fname)
        read(unit=10,fmt=*) ndesc
        read(unit=10,fmt=*) npoints
        print *,ndesc,npoints
        allocate(d(npoints))
        do i=1,npoints
            allocate(d(i)%descriptor(ndesc-1))
            read(unit=10,fmt=*) d(i)%x,d(i)%y,d(i)%a,d(i)%b,d(i)%c,d(i)%l,(d(i)%descriptor(j),j=1,ndesc-1)
        end do
        close(unit=10)
    end subroutine surf_read

    subroutine surf_file(d,fname)
        character(*) fname
        type(surfdesc), allocatable :: d(:)

        call system("/usr/local/surf/surf.ln -i "//fname//" -o /tmp/surf")
        call surf_read(d,"/tmp/surf")
    end subroutine surf_file

    subroutine surf(d,image)
        use jpegio
        type(surfdesc), allocatable :: d(:)
        real image(:,:)

        call write_jpeg2("/tmp/surf.jpg",image)
        call system("convert /tmp/surf.jpg /tmp/surf.pgm")
        call surf_file(d,"/tmp/surf.pgm")
    end subroutine surf
        
    subroutine test1
        use jpegio
        real, allocatable :: image(:,:)
        type(surfdesc), allocatable :: d(:)

        call read_jpeg2("lenna.jpg",image)
        image = 255.0 * image / maxval(image)
        call surf(d,image)
        do i=1,size(d,1)
            print *,i,d(i)%x,d(i)%y
        end do
    end subroutine test1
end module surfmod
