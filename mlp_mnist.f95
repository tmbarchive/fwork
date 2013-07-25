program mlp_mnist
  use dataio
  use mlps
  type(mlp) net
  type(autoparams) ps
  real, allocatable :: inputs(:,:),classes(:)
  real, allocatable :: test_inputs(:,:),test_classes(:)
  real errs,rate
  integer nclasses

  call read_mnist(inputs,"mnist/train-images-idx3-ubyte")
  call read_mnist(classes,"mnist/train-labels-idx1-ubyte")

  inputs = inputs / maxval(inputs)
  classes = classes + 1
  nclasses = maxval(classes)

  ps%stagesize = igetenv("stagesize",100000)
  ps%nstages = igetenv("nstages",20)
  ps%nnets = igetenv("nnets",8)
  ps%hidden_lo = fgetenv("hidden_lo",20.0)
  ps%hidden_hi = fgetenv("hidden_hi",60.0)
  ps%hidden_max = fgetenv("hidden_max",200.0)
  call mlp_autotrain_cv(ps,net,floor(classes),inputs)

  call read_mnist(test_inputs,"mnist/t10k-images-idx3-ubyte")
  call read_mnist(test_classes,"mnist/t10k-labels-idx1-ubyte")

  test_inputs = test_inputs / maxval(test_inputs)
  test_classes = test_classes + 1

  errs = mlp_error(net,floor(test_classes),test_inputs)
  rate = errs * 1.0/size(test_inputs,1)
  print *,"test","errs",errs,"total",size(test_inputs,1),"rate",rate
end program mlp_mnist
