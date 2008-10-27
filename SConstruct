# -*- Python -*-

opt = " -g -fbounds-check"
# opt = " -g -fopenmp -fbounds-check -O4"

sources = ["dataio.f95","mlpclass.f95"]

env = Environment(FORTRAN="gfortran"+opt,
                  F90="gfortran"+opt,
                  F95="gfortran"+opt,
                  LINK="gfortran"+opt)
                  
env.Program("mlpclass",sources)

