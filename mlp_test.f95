program mlp_test
  use mlpmod
  type(mlp) net
  integer, parameter :: n = 100
  real v(n,n)
  integer i,j,errs
  real err
  call mlp_init(net,n,n/2,n)
  v = 0
  forall (i=1:n) v(i,i) = 1
  do i=1,20
     print *,"***",i
     errs = 0
     err = 0
     call mlp_clear_info(net)
     do j=1,100
        v = 0
        v(1,modulo(j,n)+1) = 1
        call mlp_train(net,v,v,0.1)
     end do
     print *,i,mlp_error(net,v,v)
     call mlp_print_info(net)
  end do
end program mlp_test

