# -*- Python -*-

fc = "gfortran"
opt = " -O4"
sources = ["dataio.f95","mlpmod.f95","knnmod.f95","quicksort.f95",
           "quantizers.f95","knneval.f95","jpegio.f95","jpegio_c.c"]

env = Environment(FORTRAN=fc+opt,
                  F90=fc+opt,
                  F95=fc+opt,
                  LINK=fc+opt,
                  LIBS="-ljpeg",
                  CFLAGS="-g")
                  
env.Program("mlp_mnist",["mlp_mnist.f95"]+sources)
