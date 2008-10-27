module mlpclass

  implicit none
  
  type mlp
     integer ninput,nhidden,noutput
     real, allocatable :: b1(:), w1(:,:), b2(:), w2(:,:)
  end type mlp

  ! for now, use global variable to make wrapping easier

  type(mlp) :: net

  real :: dsigmoid_floor = 0.0

contains

  ! decay
  ! add noise
  ! normalize
  ! adjust rates

  pure real function sigmoid(x)
    real, intent(in) :: x
    sigmoid = 1.0 / (1.0+exp(-min(max(x,-200.0),200.0)))
  end function sigmoid

  pure real function dsigmoidy(y)
    real, intent(in) :: y
    dsigmoidy = max(y*(1-y),dsigmoid_floor)
  end function dsigmoidy

  pure real function dist(u,v)
    real, intent(in) :: u(:),v(:)
    dist = sqrt(sum((u-v)**2))
  end function dist

  integer function randint(n)
    integer n
    real v
    call random_number(v)
    randint = modulo(floor(n*v),n)+1
  end function randint

  subroutine mlp_init(ninput,nhidden,noutput)
    integer ninput,nhidden,noutput
    net%ninput = ninput
    net%nhidden = nhidden
    net%noutput = noutput
    allocate(net%w1(nhidden,ninput))
    allocate(net%b1(nhidden))
    allocate(net%w2(noutput,nhidden))
    allocate(net%b2(noutput))
    call random_number(net%w1)
    call random_number(net%b1)
    call random_number(net%w2)
    call random_number(net%b2)
  end subroutine mlp_init

  subroutine mlp_forward(z,x)
    real, intent(out) :: z(:)
    real, intent(in) :: x(:)
    real y(net%nhidden)
    integer i
    y = matmul(net%w1,x) + net%b1
    forall (i=1:net%nhidden) y(i) = sigmoid(y(i))
    z = matmul(net%w2,y) + net%b2
    forall (i=1:net%noutput) z(i) = sigmoid(z(i))
  end subroutine mlp_forward

  function mlp_error(actuals,xs) result(errs)
    real, intent(in) :: actuals(:,:)
    real, intent(in) :: xs(:,:)
    real :: z(size(actuals,2))
    real :: errs(2)
    integer i,row
    if (size(xs,1)/=size(actuals,1)) stop "input size mismatch"
    if (size(xs,2)/=net%ninput) stop "input vector size mismatch"
    if (size(actuals,2)/=net%noutput) stop "output vector size mismatch"
    errs = 0
    do row=1,size(xs,1)
       call mlp_forward(z,xs(row,:))
       errs(1) = errs(1) + dist(z,actuals(row,:))
       if (maxloc(z,1) /= maxloc(actuals(row,:),1)) errs(2) = errs(2) + 1
    end do
  end function mlp_error

  subroutine mlp_train(actuals,xs,eta)
    real, intent(in), target :: actuals(:,:)
    real, intent(in), target :: xs(:,:),eta
    real x(net%ninput),actual(net%noutput)
    real z(net%noutput)
    real y(net%nhidden)
    real delta2(net%noutput),delta1(net%nhidden)
    integer i,j,row

    if (size(xs,1)/=size(actuals,1)) stop "input size mismatch"
    if (size(xs,2)/=net%ninput) stop "input vector size mismatch"
    if (size(actuals,2)/=net%noutput) stop "output vector size mismatch"

    do row=1,size(xs,1)
       x = xs(row,:)
       actual = actuals(row,:)
       
       if (maxval(x)>10 .or. minval(x)<-10) stop "should normalize mlp inputs"
       ! forward propagation
       y = matmul(net%w1,x) + net%b1
       forall (i=1:net%nhidden) y(i) = sigmoid(y(i))
       z = matmul(net%w2,y) + net%b2
       forall (i=1:net%noutput) z(i) = sigmoid(z(i))

       ! backward propagation of deltas
       forall (i=1:net%noutput) delta2(i) = (z(i)-actual(i)) * dsigmoidy(z(i))
       ! delta1 = matmul(transpose(net%w2),delta2)
       delta1 = matmul(delta2,net%w2)
       forall (i=1:net%nhidden) delta1(i) = delta1(i) * dsigmoidy(y(i))
       !print *,delta2
       ! weight update
       forall (i=1:net%noutput,j=1:net%nhidden) 
          net%w2(i,j) = net%w2(i,j) - eta * delta2(i) * y(j)
       end forall
       forall (i=1:net%noutput) net%b2(i) = net%b2(i) - eta * delta2(i)
       forall (i=1:net%nhidden,j=1:net%ninput) 
          net%w1(i,j) = net%w1(i,j) - eta * delta1(i) * x(j)
       end forall
       forall (i=1:net%nhidden) net%b1(i) = net%b1(i) - eta * delta1(i)
    end do
  end subroutine mlp_train

  subroutine mlp_decay(delta)
    real delta
    integer i,j
    if (delta<0.5 .or. delta>1.0) stop "bad delta value in mlp_decay"
    forall (j=1:net%noutput,i=1:net%nhidden) 
       net%w2(i,j) = net%w2(i,j) * delta
    end forall
    forall (j=1:net%nhidden,i=1:net%ninput) 
       net%w1(i,j) = net%w1(i,j) * delta
    end forall
  end subroutine mlp_decay

  subroutine mlp_test
    integer, parameter :: n = 4
    real v(n,n)
    integer i,j,errs
    real err
    call mlp_init(n,n/2,n)
    v = 0
    forall (i=1:n) v(i,i) = 1
    do i=1,1000
       errs = 0
       err = 0
       do j=1,1000
          v = 0
          v(1,modulo(j,n)+1) = 1
          call mlp_train(v,v,0.1)
       end do
       print *,i,mlp_error(v,v)
    end do
  end subroutine mlp_test
end module

program main
  use mlpclass
  call mlp_test()
end program main
