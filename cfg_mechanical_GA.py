'''this is the code for the DHN model involved with mechanical pain network
ADVANCED VERSION WITH DR.PRESCOTT'S ADVICE
created by K. Sekiguchi (15thAUG20)
'''

from netpyne import specs
from neuron import h

cfg = specs.SimConfig()  

cfg.hParams = {'celsius': 36, 'v_init': -60 }
cfg.vrest = cfg.hParams['v_init']
cfg.duration = 5000

cfg.recordStims = False  
cfg.recordStep = 0.075# 0.025

#*--------------------------------*#
#*--- PARAMETERS FOR NETPARAMS ---*#
#*--------------------------------*#

### STIMULATION RATIO OF INPUT FIBERS ###
cfg.stim_ratios = 0.5 # 50mN -> 0.125, 100mN -> 0.25, 200mN -> 0.5, 400mN -> 1.0
if cfg.stim_ratios == 1.0: cfg.freq = '400mN'
elif cfg.stim_ratios == 0.5: cfg.freq = '200mN'
elif cfg.stim_ratios == 0.25:cfg.freq = '100mN'
else: cfg.freq = '50mN'

# SYNAPTIC WEIGHT
cfg.Ab_EX_AMPA = 0.45448           # FIXED
cfg.Ab_EX_NMDA = 0.23692           # FIXED
cfg.Ab_IN_AMPA = 0.042731          # FIXED
cfg.Ab_IN_NMDA = 0.20141           # FIXED
cfg.Ad_AMPA = 2.00e-05             # TEMPORARILY SET
cfg.Ad_NMDA = 2.00e-05             # TEMPORARILY SET
cfg.C_EX_AMPA = 0.0012924 #3.0e-03
cfg.C_EX_NMDA = 0.0017760 #3.0e-03
cfg.C_TrC_AMPA = 0.43921           # TEMPORARILY FIXED
cfg.C_DYN_AMPA = 0.49103           # TEMPORARILY FIXED
cfg.C_ISLET_AMPA = 0.32943         # TEMPORARILY FIXED
cfg.C_NK1_AMPA = 2.5108E-06 #6.0e-7
cfg.C_NK1_NMDA = 1.4946E-06 #6.0e-7
cfg.C_NK1_NK1 = 6.5293E-07 #6.0e-7

cfg.VGLUT3_PKC_AMPA = 0.11086           # FIXED
cfg.VGLUT3_PKC_NMDA = 0.10366           # FIXED
cfg.PV_GABA = 0.39221                   # FIXED
cfg.PV_GLY =  0.015361                  # FIXED
cfg.DYN_ISLET_GABA = 0.48243            # TEMPORARILY FIXED
cfg.ISLET_GABA = 0.45724                # TEMPORARILY FIXED
cfg.DYN_EX_GABA = 0.00015101 #9.0e-5
cfg.DYN_EX_GLY = 0.00014987 #9.0e-5
cfg.PKC_AMPA = 0.0011603 #8.0e-4
cfg.PKC_NMDA = 0.0023425 #8.0e-4
cfg.TrC_AMPA = 0.003
cfg.TrC_NMDA = 0.004
cfg.VGLUT3_SOM_AMPA = 0.0018624 #8.0e-4
cfg.VGLUT3_SOM_NMDA = 0.0015658 #8.0e-4
cfg.DOR_AMPA = 0.025967
cfg.DOR_NMDA = 0.014464
cfg.EX_NK1_AMPA = 0.000025598 #8.0e-7
cfg.EX_NK1_NMDA = 8.8367E-06 #8.0e-7
cfg.EX_NK1_NK1 = 0.000017163 #8.0e-7
cfg.DYN_NK1_GABA = 0.00019865
cfg.DYN_NK1_GLY =  0.00026052

cfg.recordTraces['vs'] = {'sec':'soma', 'loc':0.5,'var':'v'}

# Saving
cfg.simLabel = 'test_GA'
cfg.saveFolder = 'data'
cfg.savePickle = False
cfg.saveJson = True
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams']

# Analysis and plotting 

# script for plotting of DHN
cells = [x for x in range(442, 476, 1)]
# cfg.analysis['plotTraces'] = {'include': ['all'], 'timeRange': [0, cfg.duration], 'oneFigPer': 'trace', 'overlay': False, 'showFig' : True, 'saveFig': False} # plot recorded traces for this list of cells

cfg.analysis['plotRaster'] = {'include': ['all'], 'timeRange': [0, cfg.duration], 'saveFig': True, 'showFig': False} #'raster.png'
# cfg.analysis['plotSpikeHist'] = {'include': ['eachPop'], 'timeRange': [0,cfg.duration], 'spikeHistBin': 5, 'saveFig': True, 'showFig': False}
# cfg.analysis['plotSpikeStats'] = {'include': ['eachPop'], 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False}
# cfg.analysis['plotConn'] = {'includePre': ['all'], 'includePost': ['all'], 'feature': 'strength', 'saveFig': True, 'showFig': False}
# cfg.analysis['plotTraces'] = {'include': cells, 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False}

cfg.analysis['plot2Dnet'] = False 

# use for GA simulation
cfg.verbose = True
cfg.filename = 'sim_EXP_noise04'
cfg.printPopAvgRates = True
cfg.dt = 0.075 # 0.025

'''list of synaptic weight having all different values
cfg.synapse_ratio = 1.0
cfg.Ab_PKC_AMPA = 0.0
cfg.Ab_PKC_NMDA = 0.0
cfg.Ab_VGLUT3_AMPA = 0.0
cfg.Ab_VGLUT3_NMDA = 0.0
cfg.Ab_DOR_AMPA = 0.0
cfg.Ab_DOR_NMDA = 0.0
cfg.Ab_PV_AMPA = 0.0
cfg.Ab_PV_NMDA = 0.0
cfg.Ab_DYN_AMPA = 0.0
cfg.Ab_DYN_NMDA = 0.0
cfg.Ad_SOM_AMPA = 0.0
cfg.Ad_SOM_NMDA = 0.0
cfg.Ad_CR_AMPA = 0.0
cfg.Ad_CR_NMDA = 0.0
cfg.Ad_DOR_AMPA = 0.0
cfg.Ad_DOR_NMDA = 0.0
cfg.Ad_ISLET_GABA = 0.0
cfg.C_SOM_AMPA = 0.0
cfg.C_SOM_NMDA = 0.0
cfg.C_CR_AMPA = 0.0
cfg.C_CR_NMDA = 0.0
cfg.C_CENTRAL_AMPA = 0.0
cfg.C_ISLET_AMPA = 0.0
cfg.C_NK1_AMPA = 0.0
cfg.C_NK1_NMDA = 0.0
cfg.C_NK1_NK1 = 0.0
cfg.VGLUT3_PKC_AMPA = 0.0
cfg.VGLUT3_PKC_NMDA = 0.0
cfg.PV_PKC_GABA = 0.0
cfg.PV_PKC_GLY = 0.0
cfg.PV_DOR_GABA = 0.0
cfg.PV_DOR_GLY = 0.0
cfg.PV_CENTRAL_AMPA = 0.0
cfg.PV_CENTRAL_NMDA = 0.0
cfg.ISLET_CENTRAL_GABA = 0.0
cfg.CENTRAL_SOM_AMPA = 0.0
cfg.CENTRAL_SOM_NMDA = 0.0
cfg.PKC_SOM_AMPA = 0.0
cfg.PKC_SOM_NMDA = 0.0
cfg.VGLUT3_SOM_AMPA = 0.0
cfg.VGLUT3_SOM_NMDA = 0.0
cfg.DOR_SOM_AMPA = 0.0
cfg.DOR_SOM_NMDA = 0.0
cfg.DYN_SOM_GABA = 0.0
cfg.DYN_SOM_GLY = 0.0
cfg.CENTRAL_CR_AMPA = 0.0
cfg.CENTRAL_CR_NMDA = 0.0
cfg.PKC_CR_AMPA = 0.0
cfg.PKC_CR_NMDA = 0.0
cfg.DOR_CR_AMPA = 0.0
cfg.DOR_CR_NMDA = 0.0
cfg.DYN_CR_GABA = 0.0
cfg.DYN_CR_GLY = 0.0
cfg.SOM_NK1_AMPA = 0.0
cfg.SOM_NK1_NMDA = 0.0
cfg.SOM_NK1_NK1 = 0.0
cfg.CR_NK1_AMPA = 0.0
cfg.CR_NK1_NMDA = 0.0
cfg.CR_NK1_NK1 = 0.0
cfg.DYN_NK1_GABA = 0.0
cfg.DYN_NK1_GLY = 0.0
'''
