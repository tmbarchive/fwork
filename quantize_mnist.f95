program mlp_mnist
  use dataio
  use quantize
  type(quantparams) qs
  type(quantizer) quant
  real, allocatable :: inputs(:,:),classes(:)
  real, allocatable :: test_inputs(:,:),test_classes(:)

  call read_mnist(inputs,"mnist/train-images-idx3-ubyte")
  call read_mnist(classes,"mnist/train-labels-idx1-ubyte")

  inputs = inputs / maxval(inputs)
  classes = classes + 1
  nclasses = maxval(classes)

  call quant_train(qs,quant,inputs,10)
end program mlp_mnist
