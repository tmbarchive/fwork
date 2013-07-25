program EX1
  use jpegio
  integer PGOPEN, I
  integer w,h
  real image(100,100)
  forall(i=1:100,j=1:100) image(i,j) = i*j/10000.0
  w = 100
  h = 100
  print *,minval(image),maxval(image)
  if (PGOPEN('?') .lt. 1) stop
  call PGGRAY(image,w,h,1,w,1,h,&
    &maxval(image),minval(image),[0.0,1.0,0.0,0.0,0.0,1.0])
  call PGCLOS
  print *,"done"
end program EX1
