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

  subroutine mlp_backward(z,x,actual,eta)
    real, intent(out) :: z(:),actual(:)
    real, intent(in) :: x(:),eta
    real y(net%nhidden),delta2(net%noutput),delta1(net%nhidden)
    integer i,j

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
  end subroutine mlp_backward

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
    integer, parameter :: n = 20
    real v(n),y(n)
    integer i,j,errs
    real err
    call mlp_init(n,n/2,n)
    do i=1,1000
       errs = 0
       err = 0
       do j=1,10000
          v = 0
          v(modulo(j,n)+1) = 1
          call mlp_backward(y,v,v,1.0)
          if (maxloc(v,1) /= maxloc(y,1)) errs = errs+1
          err = err + dist(y,v)
       end do
       print *,i,errs
       print *,i,err/10000
    end do
  end subroutine mlp_test
end module

program main
  use mlpclass
  call mlp_test()
end program main
