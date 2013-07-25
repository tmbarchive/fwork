module jpegio
  implicit none
  interface
     subroutine c_read_jpeg(name,image,width,height) bind(c)
       use iso_c_binding
       integer(c_int8_t) image(*)
       integer(c_int) width,height
       character(c_char) name(*)
     end subroutine c_read_jpeg
     subroutine c_write_jpeg(name,quality,image,width,height) bind(c)
       use iso_c_binding
       integer(c_int) quality
       integer(c_int8_t) image(*)
       integer(c_int) width,height
       character(c_char) name(*)
     end subroutine c_write_jpeg
  end interface

  integer :: verbose_jpeg = 0
  integer :: jpeg_quality = 100

contains

  ! read a color JPEG image as a gray scale image

  subroutine read_jpeg2(name,image)
    character(*) name
    real, allocatable :: image(:,:)
    real, allocatable :: temp(:,:,:)
    call read_jpeg(name,temp)
    allocate(image(size(temp,2),size(temp,3)))
    image = sum(temp,1)/3.0
  end subroutine read_jpeg2

  ! write a grayscale image as a color JPEG image

  subroutine write_jpeg2(name,image)
    character(*) name
    real :: image(:,:)
    real, allocatable :: temp(:,:,:)
    logical unnormalized
    allocate(temp(3,size(image,1),size(image,2)))
    temp(1,:,:) = image
    temp(2,:,:) = image
    temp(3,:,:) = image
    call write_jpeg(name,temp)
  end subroutine write_jpeg2

  ! read a JPEG image

  subroutine read_jpeg(name,image)
    character(*) name
    real, allocatable :: image(:,:,:)
    integer(1) buffer(10000000)
    integer w,h
    call c_read_jpeg(name//char(0),buffer,w,h)
    if (verbose_jpeg>0) print *,"[read",w,h,"jpeg image]"
    allocate(image(3,w,h))
    image = reshape(buffer(1:3*w*h),[3,w,h])
    where (image<0) image = image+256
  end subroutine read_jpeg

  ! write a JPEG image

  subroutine write_jpeg(name,image)
    character(*) name
    real :: image(:,:,:)
    integer(1) buffer(10000000)
    integer w,h
    if(size(image,1)/=3) stop "images must be [3,w,h]"
    w = size(image,2)
    h = size(image,3)
    if (verbose_jpeg>0) print *,"[writing",w,h,"jpeg image]"
    buffer(1:3*w*h) = floor(reshape(image,[3*w*h]))
    call c_write_jpeg(name//char(0),jpeg_quality,buffer,w,h)
  end subroutine write_jpeg

  subroutine rgb2channels(image)
    real, allocatable :: image(:,:,:)
    real, allocatable :: temp(:,:,:)
    integer i
    allocate(temp(size(image,2),size(image,3),size(image,1)))
    forall (i=1:size(image,1)) temp(:,:,i) = image(i,:,:)
    deallocate(image)
    allocate(image(size(temp,1),size(temp,2),size(temp,3)))
    image = temp
  end subroutine rgb2channels

  subroutine channels2rgb(image)
    real, allocatable :: image(:,:,:)
    real, allocatable :: temp(:,:,:)
    integer i
    allocate(temp(size(image,3),size(image,1),size(image,2)))
    forall (i=1:size(image,3)) temp(i,:,:) = image(:,:,i)
    deallocate(image)
    allocate(image(size(temp,1),size(temp,2),size(temp,3)))
    image = temp
  end subroutine channels2rgb
  
  subroutine test_ll
    integer(1) image(100000000)
    integer w,h
    call c_read_jpeg("test.jpg"//char(0),image,w,h)
    image(1:100:3) = 0
    image(2:100:3) = -127
    image(3:100:3) = 0
    call c_write_jpeg("test2.jpg"//char(0),90,image,w,h)
    print *,w,h
  end subroutine test_ll
  
  subroutine test_hl
    real, allocatable :: temp(:,:,:),image(:,:,:)
    call read_jpeg("test.jpg",temp)
    call rgb2channels(temp)
    temp(1:10,1:17,2) = 255
    call channels2rgb(temp)
    call write_jpeg("test-x.jpg",temp)
  end subroutine test_hl
end module jpegio

