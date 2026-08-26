c Common blocks for NLO or aMC@NLO reweighting
c Used only for interfacing with APPLgrid/PineAPPL
c Counterevents use icountevts=0,1,2 for the soft, collinear, and
c soft-collinear contexts, respectively.

c The following is introduced since this file may be included in codes
c that need not know nexternal. Make sure that maxparticles>=nexternal
      integer maxparticles
      parameter (maxparticles=20)
