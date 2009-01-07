module dataio

    implicit none

    interface read_mnist
        module procedure read_mnist1, read_mnist2, read_mnist3,&
            &read_mnist1i, read_mnist2i, read_mnist3i
    end interface

    interface display
        module procedure display_pgm, display_ppm, display_1d
    end interface

contains

    !----------------------------------------------------------------
    ! reading/writing number tables
    !----------------------------------------------------------------

    function countfields(line) result(count)
        ! count the number of non-blank fields in the line
        character(*) line
        integer count,i,n,last
        n = len(line)
        i = 1
        ! skip initial blanks
        do while (i<=n .and. line(i:i)==" ")
            i=i+1;
        end do
        count = 0
        ! iterate through the fields
        do while(i<n)
            count = count+1
            do while (i<n .and. line(i:i)/=" ")
                i=i+1
            end do
            do while (i<n .and. line(i:i)==" ")
                i=i+1
            end do
        end do
    end function countfields

    ! read a text file into an allocatable array
    ! the text file contains blank-separated numerical fields

    subroutine read_table(a,fname)
        ! Fortran is still sadly deficient in memory
        ! management and I/O... unit=, no realloc, ...
        real,allocatable,intent(out) :: a(:,:)
        character(*) fname
        character(1000000) line
        integer count,rows,cols,row,col

        ! first pass: determine size
        open(unit=9,file=fname)
        rows = 0
        cols = 0
        do
            read(unit=9,fmt='(A)',end=900) line
            rows = rows+1
            ! get number of columns from first non-zero line
            if (cols==0) cols = countfields(line)
        end do
900     continue

        ! second pass: allocate output, then read the data into it
        rewind(unit=9)
        allocate(a(rows,cols))
        a(:,:) = 0
        row=0
        do
            row=row+1
            read(unit=9,fmt='(A)',end=901) line
            read(line,fmt=*) (a(row,col),col=1,cols)
        end do
901     continue
        close(unit=9)
    end subroutine read_table

    subroutine write_table(fname,a)
        real a(:,:),lo,hi
        character(*) fname
        integer i,j
        open(10,file=trim(fname))
        do i=1,size(a,1)
            write(10,fmt=*) (a(i,j),j=1,size(a,2))
        end do
        close(10)
    end subroutine write_table

    !----------------------------------------------------------------
    ! reading/writing MNIST data
    !----------------------------------------------------------------

    integer function readbyte(u)
        integer u
        integer(1) byte
        read(unit=u)byte
        readbyte = byte
        if (readbyte<0) readbyte=256+readbyte
    end function readbyte

    integer function readint32(u)
        integer u,val,i
        val = 0
        do i=1,4
            val = val*256+readbyte(u)
        end do
        readint32 = val
    end function readint32

    subroutine read_mnist1(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:)
        real, allocatable :: a(:)
        integer nitems,magic,i
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic/=2049) stop "bad magic number"
        nitems = readint32(10)
        print *,trim(fname),nitems
        allocate(temp(nitems))
        read(10)temp
        allocate(a(nitems))
        forall (i=1:nitems) a(i) = temp(i)
        where (a<0) a = a+256
    end subroutine read_mnist1

    subroutine read_mnist2(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:,:),temp3(:,:,:)
        real, allocatable :: a(:,:)
        integer magic,rows,cols,i,j,k,nitems
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic==2050) then
            rows = readint32(10)
            cols = readint32(10)
            print *,trim(fname),rows,cols
            allocate(temp(cols,rows))
            read(10)temp
            allocate(a(rows,cols))
            forall (i=1:rows,j=1:cols) a(i,j) = temp(j,i)
            where (a<0) a = a+256
        else if(magic==2051) then
            ! common special case: read rank 3 data on disk into rank 2 array
            nitems = readint32(10)
            rows = readint32(10)
            cols = readint32(10)
            print *,trim(fname),nitems,rows,cols
            allocate(temp3(cols,rows,nitems))
            read(10)temp3
            allocate(a(nitems,rows*cols))
            do i=1,nitems
                do j=1,rows
                    do k=1,cols
                        a(i,(j-1)*cols+k) = temp3(k,j,i)
                    end do
                end do
            end do
            where (a<0) a = a+256
        else
            stop "bad magic number"
        end if
    end subroutine read_mnist2

    subroutine read_mnist3(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:,:,:)
        real, allocatable :: a(:,:,:)
        integer nitems,magic,rows,cols,i,j,k
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic/=2051) stop "bad magic number in read_mnist3"
        nitems = readint32(10)
        rows = readint32(10)
        cols = readint32(10)
        print *,trim(fname),nitems,rows,cols
        allocate(temp(cols,rows,nitems))
        read(10)temp
        allocate(a(nitems,rows,cols))
        forall (i=1:nitems,j=1:rows,k=1:cols) a(i,j,k) = temp(k,j,i)
        where (a<0) a = a+256
    end subroutine read_mnist3

    subroutine read_mnist1i(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:)
        integer, allocatable :: a(:)
        integer nitems,magic,i
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic/=2049) stop "bad magic number"
        nitems = readint32(10)
        print *,trim(fname),nitems
        allocate(temp(nitems))
        read(10)temp
        allocate(a(nitems))
        forall (i=1:nitems) a(i) = temp(i)
        where (a<0) a = a+256
    end subroutine read_mnist1i

    subroutine read_mnist2i(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:,:)
        integer, allocatable :: a(:,:)
        integer magic,rows,cols,i,j
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic==2050) then
            rows = readint32(10)
            cols = readint32(10)
        else if(magic==2051) then
            rows = readint32(10)
            cols = readint32(10)*readint32(10)
        else
            stop "bad magic number"
        end if
        print *,trim(fname),rows,cols
        allocate(temp(rows,cols))
        read(10)temp
        allocate(a(rows,cols))
        forall (i=1:rows,j=1:cols) a(i,j) = temp(i,j)
        where (a<0) a = a+256
    end subroutine read_mnist2i

    subroutine read_mnist3i(a,fname)
        character(*) fname
        integer(1), allocatable :: temp(:,:,:)
        integer, allocatable :: a(:,:,:)
        integer nitems,magic,rows,cols,i,j,k
        open(unit=10,file=trim(fname),status='old',access='stream')
        magic = readint32(10)
        if (magic/=2051) stop "bad magic number in read_mnist3"
        nitems = readint32(10)
        rows = readint32(10)
        cols = readint32(10)
        print *,trim(fname),nitems,rows,cols
        allocate(temp(nitems,rows,cols))
        read(10)temp
        allocate(a(nitems,rows,cols))
        forall (i=1:nitems,j=1:rows,k=1:cols) a(i,j,k) = temp(i,j,k)
        where (a<0) a = a+256
    end subroutine read_mnist3i

    !----------------------------------------------------------------
    ! reading/writing NetPBM image data
    !----------------------------------------------------------------

    subroutine display_1d(a,w,h)
        integer w,h,i,j
        real a(:)
        real f(w,h)
        do i=1,w
            do j=1,h
                f(i,j) = a((i-1)*h+j)
            end do
        end do
        call display_pgm(f)
    end subroutine display_1d

    subroutine display_pgm(a)
        real a(:,:)
        call system("rm -f /tmp/__temp__.pgm")
        call write_pgm("/tmp/__temp__.pgm",a)
        call system("display /tmp/__temp__.pgm")
    end subroutine display_pgm

    subroutine display_ppm(a)
        real a(:,:,:)
        call system("rm -f /tmp/__temp__.pgm")
        call write_ppm("/tmp/__temp__.ppm",a)
        call system("display /tmp/__temp__.ppm")
    end subroutine display_ppm

    subroutine write_pgm(fname,a)
        character(*) fname
        integer i,j,byte
        real a(:,:),lo,hi
        lo = minval(a)
        hi = maxval(a)
        open(10,file=trim(fname))
        write(10,fmt="(a/)",advance="no") "P2"
        write(10,fmt="(i0,a,i0,a,i0/)",advance="no") size(a,1)," ",size(a,2)," ",255
        do j=1,size(a,2)
            do i=1,size(a,1)
                byte = max(0,min(255,floor(256 * (a(i,j)-lo) / (hi-lo))))
                write(10,fmt="(a)",advance="no") " "
                write(10,fmt="(i0)",advance="no") byte
            end do
            write(10,fmt="(a/)",advance="no") ""
        end do
        close(1)
    end subroutine write_pgm

    subroutine write_ppm(fname,a)
        character(*) fname
        integer i,j,k,byte
        real a(:,:,:),lo,hi
        lo = minval(a)
        hi = maxval(a)
        open(10,file=trim(fname))
        write(10,fmt="(a/)",advance="no") "P3"
        write(10,fmt="(i0,a,i0,a,i0/)",advance="no") size(a,1)," ",size(a,2)," ",255
        do j=1,size(a,2)
            do i=1,size(a,1)
                do k=1,size(a,3)
                    byte = max(0,min(255,floor(256 * (a(i,j,k)-lo) / (hi-lo))))
                    write(10,fmt="(a)",advance="no") " "
                    write(10,fmt="(i0)",advance="no") byte
                end do
            end do
            write(10,fmt="(a/)",advance="no") ""
        end do
        close(1)
    end subroutine write_ppm

    !----------------------------------------------------------------
    ! simple test cases
    !----------------------------------------------------------------

    subroutine test_mnist
        real, allocatable :: a(:,:,:), l(:)
        call read_mnist(a,"mnist/train-images-idx3-ubyte")
        call read_mnist(l,"mnist/train-labels-idx1-ubyte")
    end subroutine test_mnist

    subroutine test_table
        real, allocatable :: a(:,:)
        character(100) fname
        call getarg(1,fname)
        print *,fname
        call read_table(a,fname)
        print *,size(a,1),size(a,2)
        print *,a
    end subroutine test_table

end module dataio
