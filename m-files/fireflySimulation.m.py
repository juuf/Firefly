# Automatically translated using m2cpp 2.0 on 2019-04-12 16:39:20
#
# encoding: utf-8
#
# Supplement file
#
# Valid inputs:
#
# uword   int     float   double cx_double
# uvec    ivec    fvec    vec    cx_vec
# urowvec irowvec frowvec rowvec cx_rowvec
# umat    imat    fmat    mat    cx_mat
# ucube   icube   fcube   cube   cx_cube
#
# char    string  struct  structs func_lambda

functions = {
  "fireflySimulation" : {
                       "Alpha" : "",
                        "Beta" : "",
           "DoublingThreshold" : "",
                           "E" : "", # mat
                "FilterLength" : "",
                  "PlotLength" : "",
            "SimulationLength" : "",
                   "animSpeed" : "",
                     "animate" : "",
                "defaultAlpha" : "", # double
            "defaultAnimSpeed" : "", # double
                 "defaultBeta" : "", # double
    "defaultDoublingThreshold" : "", # int
         "defaultFilterLength" : "", # int
           "defaultFilterType" : "", # string
                "defaultNodes" : "", # int
           "defaultPlotLength" : "", # int
           "defaultSampleRate" : "", # int
                 "defaultSeed" : "", # int
     "defaultSimulationLength" : "", # int
                 "defaulttref" : "", # int
                      "domega" : "", # mat
                           "f" : "", # rowvec
                  "filterType" : "",
                        "fire" : "", # int
                           "i" : "", # int
                          "in" : "", # int
                "interactions" : "",
                           "j" : "", # int
                           "k" : "", # int
                          "k2" : "",
                          "k3" : "",
                          "k4" : "",
                          "l1" : "",
                          "l3" : "",
                          "l4" : "",
                           "n" : "",
                       "omega" : "", # mat
                    "omegaOut" : "",
                      "omegas" : "", # mat
                           "p" : "",
                          "p1" : "",
                          "p2" : "",
                          "pf" : "", # int
                    "phaseOut" : "",
                         "phi" : "", # mat
                      "phiref" : "",
             "plotSpectrogram" : "",
                           "q" : "", # int
                           "r" : "",
                  "randomSeed" : "",
                "sortedOmegas" : "",
                          "sr" : "",
                         "tmp" : "",
                        "tref" : "",
                    "varargin" : "",
  },
}
structs = {
  "p" : {
    "KeepUnmatched" : "", # int
          "Results" : "",
  },
}
includes = [
  '#include "mconvert.h"',
  '#include <iostream>',
  '#include <armadillo>',
  'using namespace arma ;',
]