# set MaxEvents 10000
set RandomSeed 123
# set RandomSeed 123


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
  RHadronFilter
  RHadronMerger

  PileUpMerger
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  TrackSmearing

  TimeSmearing

  VertexFinderDAClusterizerZT
  HighMassVertexRecover

  TreeWriter
}

#################################
# Get stable R-Hadrons
#################################

module PdgCodeFilter RHadronFilter {
  set InputArray Delphes/allParticles
  set OutputArray stableRHadrons
  set Invert false

  set RequireStatus true
  set Status 104
}

#################################
# Merge R-Hadrons with other sable
#################################

module Merger RHadronMerger {
  add InputArray RHadronFilter/stableRHadrons
  add InputArray Delphes/stableParticles
  set OutputArray stableParticles
}


###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  # set InputArray Delphes/stableParticles
  set InputArray RHadronMerger/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  set Verbose 0

  # pre-generated minbias input file
  set PileUpFile /afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/_hepmc/MinBias_pp14TeV_100k.pileup

  # average expected pile up
  set MeanPileUp 140

  # 0-poisson, 1-uniform, 2-delta
  set PileUpDistribution 2

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.25

  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)

  #set VertexDistributionFormula {exp(-(t^2/(2*(0.063/2.99792458E8*exp(-(z^2/(2*(0.063)^2))))^2)))}
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}

  # taking 5.3 cm x 160 ps

  #set VertexDistributionFormula { (abs(t) <= 160e-12) * (abs(z) <= 0.053) * (1.00) +
  #                                (abs(t) >  160e-12) * (abs(z) <= 0.053) * (0.00) +
  #                   	          (abs(t) <= 160e-12) * (abs(z) > 0.053)  * (0.00) +
  # 				                  (abs(t) >  160e-12) * (abs(z) > 0.053)  * (0.00)}

}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.0

  # magnetic field
  set Bz 3.8
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.2)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.2   && pt <= 1.0)   * (0.70) +
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.2   && pt <= 1.0)   * (0.60) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.85) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.73) +
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.95) +
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) +
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e3) * (0.99) +
                                           (abs(eta) <= 1.5) * (pt > 1.0e3 )               * (0.99 * exp(0.5 - pt*5.0e-4)) +

                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e3) * (0.98) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e3)                * (0.98 * exp(0.5 - pt*5.0e-4)) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}


########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons ,

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackMomentumResolution.tcl
}


########################################
# Momentum resolution for electrons
########################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackMomentumResolution.tcl
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons

  # up to |eta| < 2.8 take measurement from tracking + muon chambers
  # for |eta| > 2.8 and pT < 5.0 take measurement from tracking alone taken from
  # http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source muonMomentumResolution.tcl
}


##############
# Track merger
##############

module Merger TrackMerger {
  # add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

########################################
#   Smear tracks
########################################

module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks

  set OutputArray tracks
  set ApplyToPileUp true

  set ZOuterResolution .7
  set Bz 3.8

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackResolutionCMS.tcl
}

########################################
#   Time Smearing
########################################

module TimeSmearing TimeSmearing {
  set InputArray TrackSmearing/tracks
  set OutputArray tracks

  # assume 30 ps resolution for now
  set TimeResolution 3.0E-11
}

##################################
# Primary vertex clustering
##################################

module VertexFinderDAClusterizerZT VertexFinderDAClusterizerZT {
  set InputArray TimeSmearing/tracks

  set TrackOutputArray tracks
  set VertexOutputArray vertices

  set Verbose 0

  # [mm]
  set DzCutOff 40
  # in d0/sigma_d0
  set D0CutOff 1
}


######################################
# Heavy(slow) particles vertex recover
######################################

module HighMassVertexRecover HighMassVertexRecover {
  set TrackInputArray VertexFinderDAClusterizerZT/tracks
  set VertexInputArray VertexFinderDAClusterizerZT/vertices

  set TrackOutputArray tracks
  set VertexOutputArray vertices

  set Verbose 0

  # add MassHypot 0.13957
}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  # add Branch InputArray BranchName BranchClass

  add Branch PileUpMerger/stableParticles Particle GenParticle
  add Branch HighMassVertexRecover/tracks Track Track

  add Branch HighMassVertexRecover/vertices Vertex4D Vertex

  add Branch PileUpMerger/vertices GenVertex Vertex


}
