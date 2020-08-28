# This is the code to create SOM positive neuron in lamina II having initial burst firing
# created by K. Sekiguchi (30th Apr 2020)

from netpyne import specs, sim
from neuron import h

try: 
    from __main__ import cfg
except:
    from cfg_EXinitial import cfg

netParams = specs.NetParams()

###################################################################################################################################
#   Create Populations
###################################################################################################################################

netParams.popParams['EX'] = {'cellType': 'EX' , 'numCells':  4, 'cellModel': '_EX' }

cellRule  = netParams.importCellParams(label    = 'EXrule' , conds = {'cellType': 'EX'  ,'cellModel': '_EX'}, 
                                       fileName = 'EXinitialburst.tem', cellName = 'EXinitial')

## SOMA ##
cellRule['secs']['soma']['geom']['Ra'] = cfg.Ra
cellRule['secs']['soma']['geom']['cm'] = cfg.cm
cellRule['secs']['soma']['mechs']   = {'B_Na' : {'gnabar' : cfg.gnabar_B_Na_soma},
                                       'B_A'  : {'gkbar'  : 0.0},
                                       'B_DR' : {'gkbar'  : cfg.gkbar_B_DR_soma },
                                       'KDR'  : {'gkbar'  : 0.0}, 
                                       'KDRI' : {'gkbar'  : cfg.gkbar_KDRI_soma },
                                       'CaIntraCellDyn'   : {'depth': 0.1, 'cai_tau': 1.0, 'cai_inf': 5e-05},                              
                                       'iKCa' : {'gbar'   : 0.002,     'gk': 0  },
                                       'pas'  : {'g'      : cfg.g_pas, 'e': cfg.e_pas}}

## DEND ##
cellRule['secs']['hillock']['geom']['Ra'] = cfg.Ra
cellRule['secs']['hillock']['geom']['cm'] = cfg.cm

cellRule['secs']['dend']['mechs'] = {'SS'  : {'gnabar' : 0.0}, 
                                     'B_DR': {'gkbar'  : 0.0}, 
                                     'KDR' : {'gkbar'  : 0.0},
                                     'KDRI': {'gkbar'  : cfg.gkbar_KDRI_dend},
                                     'CaIntraCellDyn'  : {'depth': 0.1, 'cai_tau': 1.0, 'cai_inf': 5e-05},                              
                                     'iKCa': {'gbar'   : 0.002, 'gk': 0.0},
                                     'pas' : {'g'      : cfg.g_pas, 'e': cfg.e_pas}}

## HILLOCK ##
cellRule['secs']['dend']['geom']['Ra'] = cfg.Ra
cellRule['secs']['dend']['geom']['cm'] = cfg.cm

cellRule['secs']['hillock']['mechs'] = {'B_Na': {'gnabar' : cfg.gnabar_B_Na_hill},
                                        'B_A' : {'gkbar'  : cfg.gkbar_B_A_hill  }, 
                                        'B_DR': {'gkbar'  : cfg.gkbar_B_DR_hill }, 
                                        'KDR' : {'gkbar'  : cfg.gkbar_KDR_hill  },
                                        'KDRI': {'gkbar'  : cfg.gkbar_KDRI_hill },
                                        'pas' : {'g'      : cfg.g_pas, 'e': cfg.e_pas}}

## AXON ##
cellRule['secs']['axon']['geom']['Ra'] = cfg.Ra
cellRule['secs']['axon']['geom']['cm'] = cfg.cm

cellRule['secs']['axon']['mechs'] = {'B_Na' : {'gnabar' : 0.0},
                                     'B_A'  : {'gkbar'  : 0.0}, 
                                     'B_DR' : {'gkbar'  : 0.0}, 
                                     'KDR'  : {'gkbar'  : 0.0},
                                     'KDRI' : {'gkbar'  : 0.0},
                                     'pas'  : {'g'      : 0.0, 'e': cfg.e_pas}}

cellRule['secs']['soma']['threshold'] = -10
netParams.cellParams['EXRule'] = cellRule


###################################################################################################################################
#   Stimulation from fiber
###################################################################################################################################

# IClamp for Melnick 2004 condition
# netParams.stimSourceParams['stim'] = {'type': 'IClamp', 'delay': 50, 'dur': 500, 'amp': 0.06} 
# netParams.stimTargetParams['stim'] = {'source': 'stim', 'conds': {'pop': 'EX'}, 'sec': 'soma', 'loc': 0.55}

# IClamp for SOM paper
amp = [ 0.065, 0.050, 0.035, 0.020 ]
for i, key in enumerate(amp):
    netParams.stimSourceParams[i] = {'type': 'IClamp', 'delay': 50, 'dur': 1150, 'amp': key} 
    netParams.stimTargetParams[i] = {'source': i, 'conds': {'pop': 'EX', 'cellList':[ i ]}, 'sec': 'soma', 'loc': 0.5}

# netParams.defaultThreshold = -30

h.alpha_shift_B_Na = cfg.alpha_shift_B_Na
h.beta_shift_B_Na  = cfg.beta_shift_B_Na