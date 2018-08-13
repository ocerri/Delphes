# set MaxEvents 10000
set RandomSeed 123


#
#  Phase II - Pile-Up
#
#  Main authors: Olmo Cerri (CalTech)
#
#  Released on:
#
#  Version: v02 beta - test TrackSmearing and vertexing
#
#
#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  HardProcessFilter
  HardProcessMerger

  PileUpMerger

  TreeWriter
}

#################################
# Get stable R-Hadrons
#################################

module PdgCodeFilter HardProcessFilter {
  set InputArray Delphes/allParticles
  set OutputArray hardProcessParticles
  set Invert false

  set RequireStatus true
  set Status 3
}

#################################
# Merge hard scattering with other sable
#################################

module Merger HardProcessMerger {
  add InputArray HardProcessFilter/hardProcessParticles
  add InputArray Delphes/stableParticles
  set OutputArray stableParticles
}


###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  # set InputArray Delphes/stableParticles
  set InputArray HardProcessMerger/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  set Verbose 0

  # pre-generated minbias input file
  set PileUpFile /Users/olmo/cernbox/PID_timing_studies/_hepmc/MinBias_pp14TeV_100k.pileup
  # set PileUpFile /eos/cms/store/group/upgrade/delphes/PhaseII/MinBias_100k.pileup

  # average expected pile up
  set MeanPileUp 80

  # 0-poisson, 1-uniform, 2-delta
  set PileUpDistribution 0

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.25

  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}
  # taking 5.3 cm x 160 ps

}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  # add Branch InputArray BranchName BranchClass
  add Branch PileUpMerger/stableParticles GenParticles GenParticle

}
