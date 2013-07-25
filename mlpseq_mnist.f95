program mlp_mnist
  use dataio
  use mlps
  type(mlp) net
  type(autoparams) ps
  real, allocatable :: inputs(:,:),classes(:)
  real, allocatable :: test_inputs(:,:),test_classes(:)
  integer nclasses

  call read_mnist(inputs,"mnist/train-images-idx3-ubyte")
  call read_mnist(classes,"mnist/train-labels-idx1-ubyte")

  inputs = inputs / maxval(inputs)
  classes = classes + 1
  nclasses = maxval(classes)

  call read_mnist(test_inputs,"mnist/t10k-images-idx3-ubyte")
  call read_mnist(test_classes,"mnist/t10k-labels-idx1-ubyte")

  test_inputs = test_inputs / maxval(test_inputs)
  test_classes = test_classes + 1

  call mlp_init(net,size(inputs,2),20,nclasses)
  call mlp_autotrain(ps,net,floor(classes),inputs,floor(test_classes),test_inputs)
end program mlp_mnist
