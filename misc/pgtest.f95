program EX1
  integer PGOPEN, I
  real XS(9), YS(9), XR(101), YR(101)
  do I=1,101
     XR(I) = 0.1*(I-1)
     YR(I) = XR(I)**2*exp(-XR(I))
  end do
  do I=1,9
     XS(I) = I
     YS(I) = XS(I)**2*exp(-XS(I))
  end do
  if (PGOPEN('?') .lt. 1) stop
  call PGENV(0., 10., 0., 0.65,  0,  0)
  call PGLAB('x', 'y', 'PGPLOT Graph: y = x\u2\dexp(-x)')
  call PGLINE(101, XR, YR)
  call PGPT(9, XS, YS, 18)
  call PGCLOS
end program EX1
