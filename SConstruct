# -*- Python -*-

opt = " -g"
opt += " -I/usr/lib/fortran/modules/plplot"
fc = "gfortran -std=f2003"

if 1:
    opt += " -O3"
    opt += " -fopenmp"
else:
    fc = "g95 -O"
    opt += " -fbounds-check"
#    opt += " -fno-inline"
#    opt += " -fbacktrace"
#    opt += " -finit-real=nan -finit-integer=-999999999"

env = Environment(FORTRAN=fc+opt,
                  F90=fc+opt,
                  F95=fc+opt,
                  LINK=fc+opt,
                  CFLAGS="-g")

# basic algorithms
sources = ["dataio.f95","quicksort.f95","reallocs.f95","utils.f95"]

# interfaces to other packages
sources += ["jpegio.f95","jpegio_c.c","surfmod.f95"]

# pattern recognition algorithms
sources += ["cedges.f95","quantize.f95","nnbr.f95","fextract.f95"]

# objects / classes            
sources += ["mlps.f95","knncs.f95","uniclasses.f95"]

env.Append(LIBS=["-ljpeg"])
# env.Append(LIBS=["-lefence"])

if 0:
    sources += "plutil.f95",
    env.Append(LIBS=["-lplplotd","-lplplotf95d"])

env.Program("uni_mnist",["uni_mnist.f95"]+sources)
env.Program("quantize_mnist",["quantize_mnist.f95"]+sources)
env.Program("mlp_mnist",["mlp_mnist.f95"]+sources)
env.Program("mnist_nnbr",["mnist_nnbr.f95"]+sources)

# env.Library("fclassifier",["cinterfaces.f95"]+sources)
# env.Program("uni_mnist",["uni_mnist.f95"]+sources)
# env.Program("display_cedges",["display-cedges.f95"]+sources)
# env.Program("test_cedges",["test-cedges.f95"]+sources)
# env.Program("rast-match",["rast-match.f95"]+sources)
# env.Program("test-rast",["test-rast.f95"]+sources)
