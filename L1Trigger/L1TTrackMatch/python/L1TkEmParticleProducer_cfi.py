import FWCore.ParameterSet.Config as cms

L1TkPhotons = cms.EDProducer("L1TkEmParticleProducer",
        label = cms.string("EG"), 	# labels the collection of L1TkEmParticleProducer that is produced
					# (not really needed actually)
        L1EGammaInputTag = cms.InputTag("simCaloStage2Digis",""),
                                                # When the standard sequences are used :
                                                #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
                                                #   - for the "old stage-2" algo (2x2 clustering), use 
                                                #     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
        ETmin = cms.double( -1 ),               # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
        RelativeIsolation = cms.bool( True ),   # default = True. The isolation variable is relative if True,
                                                # else absolute.
        IsoCut = cms.double( 0.23 ),             # Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
                                                # When RelativeIsolation = False, IsoCut is in GeV.
           # Determination of the isolation w.r.t. L1Tracks :
     	L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
        ZMAX = cms.double( 25. ),       # in cm
        CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 2. ),    # in GeV
        DRmin = cms.double( 0.07),
        DRmax = cms.double( 0.30 ),
        PrimaryVtxConstrain = cms.bool( False ),  # default = False
						  # if set to True, the L1TkPrimaryVertex is used to constrain
						  # the tracks entering in the calculation of the isolatiom.
        #DeltaZConstrain = cms.bool( False ),  # default = False
                                                  # if set to True, constrain to the z of the leading
						  # track within DR < DRmax
        DeltaZMax = cms.double( 999. ),    # in cm. Used only when PrimaryVtxConstrain = True
        L1VertexInputTag = cms.InputTag("NotUsed"),     # Used only when PrimaryVtxConstrain = True
)


L1TkPhotonsTightIsol = L1TkPhotons.clone()
L1TkPhotonsTightIsol.IsoCut = cms.double( 0.10)

#### Additional collections that right now only the menu team is using - to be renamed/redefined by the EGamma group
# The important change is the EG seed -> PhaseII instead of PhaseI

L1TkPhotonsCrystal=L1TkPhotons.clone()
L1TkPhotonsCrystal.L1EGammaInputTag = cms.InputTag("L1EGammaClusterEmuProducer", "L1EGammaCollectionBXVEmulator")
L1TkPhotonsCrystal.IsoCut = cms.double(-0.1)


L1TkPhotonsHGC=L1TkPhotons.clone()
L1TkPhotonsHGC.L1EGammaInputTag = cms.InputTag("l1EGammaEEProducer","L1EGammaCollectionBXVWithCuts")
L1TkPhotonsHGC.IsoCut = cms.double(-0.1)


