'''this is the code for the DHN model involved with mechanical pain network
ADVANCED VERSION WITH DR.PRESCOTT'S ADVICE
THIS IS 2ND REVISION for GENETIC ALGORITHM
created by K. Sekiguchi (15thAUG20)
<neuron number>
1/10 scale 
Ab_SAI: 0-9 (10),   Ab_SAII: 10-19 (10),  Ad: 20-39 (20),   C_PEP: 40-119 (80), C_NP: 120-199 (80)
PKC: 200-259 (60),  VGLUT3: 260-269 (10), PV: 270-299 (30), DOR: 300-359 (60),  TrC: 360-379(20)
DYN: 380-499 (120), SOM: 500-529 (30),    CR: 530-569 (40), ISLET: 570-599 (30),NK1: 600-619(20)
1/20 scale 
Ab_SAI: 0-4 (5),   Ab_SAII: 5-9 (5),  Ad: 10-19 (10),   C_PEP: 20-59 (40), C_NP: 60-99 (40)
PKC: 100-129 (30),  VGLUT3: 130-134 (5), PV: 135-149 (15), DOR: 150-179 (30),  TrC: 180-189(10)
DYN: 190-249 (60), SOM: 250-264 (15),    CR: 265-284 (20), ISLET: 285-299 (15),NK1: 300-309(10)
'''

from netpyne import specs, sim
from neuron import h
import numpy as np
import cells
from spkt_gen import *
import json
#from cfg_mechanical_GA import cfg

try: 
    from __main__ import cfg
except:
    from cfg_mechanical_GA import cfg

netParams = specs.NetParams()

#------------------------------------------------------------------------------
# Population
#------------------------------------------------------------------------------

# ### INPUT FREQUENCY FROM AFFERENT FIBERS WITH RANDOM FACTOR ###
# rate_SAI, t = rate_SAI()
# rate_SAII, t = rate_SAII()
# rate_Ad, t = rate_Ad()
# rate_C, t = rate_C()

# spkt_SAI, spkt_SAII, spkt_Ad, spkt_C = [], [], [], []

# for i in range(5):
#     seeds = np.random.randint(50)
#     spkt_SAI_ = inh_poisson_generator(rate_SAI, t, 5000, seed=seeds)
#     spkt_SAII_ = inh_poisson_generator(rate_SAII, t, 5000, seed=seeds)
#     spkt_SAI.append(spkt_SAI_)
#     spkt_SAII.append(spkt_SAII_)

# for i in range(10):
#     seeds = np.random.randint(100)
#     spkt_Ad_ = inh_poisson_generator(rate_Ad, t, 5000, seed=seeds)
#     spkt_Ad.append(spkt_Ad_)

# for i in range(40):
#     seeds = np.random.randint(400)
#     spkt_C_ = inh_poisson_generator(rate_C, t, 5000, seed=seeds)
#     spkt_C.append(spkt_C_)

# with open('./spkt_gen/spkt_SAI_%s.json' %(cfg.freq), 'w') as SAI: json.dump(spkt_SAI, SAI)
# with open('./spkt_gen/spkt_SAII_%s.json' %(cfg.freq), 'w') as SAII: json.dump(spkt_SAII, SAII)
# with open('./spkt_gen/spkt_Ad_%s.json' %(cfg.freq), 'w') as Ad: json.dump(spkt_Ad, Ad)
# with open('./spkt_gen/spkt_C_%s.json' %(cfg.freq), 'w') as  C: json.dump(spkt_C, C)

### INPUT FREQUENCY FROM AFFERENT FIBERS WITH FIXED FACTOR ###
# 50mN #
with open('spkt_SAI_200mN.json', 'rb') as spkt_SAI01: spkt_SAI01 = json.load(spkt_SAI01)
with open('spkt_SAII_200mN.json', 'rb') as spkt_SAII01: spkt_SAII01 = json.load(spkt_SAII01)
with open('spkt_Ad_200mN.json', 'rb') as spkt_Ad01: spkt_Ad01 = json.load(spkt_Ad01)
with open('spkt_C_200mN.json', 'rb') as spkt_C01: spkt_C01 = json.load(spkt_C01)

# 200mN #
with open('spkt_SAI_200mN.json', 'rb') as spkt_SAI02: spkt_SAI02 = json.load(spkt_SAI02)
with open('spkt_SAII_200mN.json', 'rb') as spkt_SAII02: spkt_SAII02 = json.load(spkt_SAII02)
with open('spkt_Ad_200mN.json', 'rb') as spkt_Ad02: spkt_Ad02 = json.load(spkt_Ad02)
with open('spkt_C_200mN.json', 'rb') as spkt_C02: spkt_C02 = json.load(spkt_C02)

cells.PKCRule['conds'] = {'cellType': 'PKC'}
cells.PKC02Rule['conds'] = {'cellType': 'PKC02'}
cells.INcellRule['conds'] = {'cellType': 'IN'}
cells.EXdelayedRule['conds'] = {'cellType': 'EXdl'}
cells.CRRule['conds'] = {'cellType': 'CR'}
cells.SOMRule['conds'] = {'cellType': 'SOM'}
cells.EXinitialRule['conds'] = {'cellType': 'EXib'}
cells.PROcellRule['conds'] = {'cellType': 'PRO'}
cells.PROcellRule['secs']['soma']['threshold'] = 0

netParams.cellParams['PKC'] = cells.PKCRule
netParams.cellParams['PKC02'] = cells.PKC02Rule
netParams.cellParams['VGLUT3Rule'] = cells.EXdelayedRule
netParams.cellParams['PVRule'] = cells.INcellRule
netParams.cellParams['DORRule'] = cells.EXdelayedRule
netParams.cellParams['TrCRule'] = cells.EXinitialRule
netParams.cellParams['DYNRule'] = cells.INcellRule
netParams.cellParams['SOMRule'] = cells.SOMRule
netParams.cellParams['CRRule'] = cells.CRRule
netParams.cellParams['ISLETRule'] = cells.INcellRule
netParams.cellParams['NK1Rule'] = cells.PROcellRule

#--** NETWORK 01 **--#
## PRIMARY AFFERENTS 
netParams.popParams['Ab_SAI01'] = {'cellModel': 'VecStim', 'numCells': 5, 'spkTimes': spkt_SAI01}  # input from Ab_slow adapting type I
netParams.popParams['Ab_SAII01'] = {'cellModel': 'VecStim', 'numCells': 5, 'spkTimes': spkt_SAII01}  # input from Ab_slow adapting type I
netParams.popParams['Ad01'] = {'cellModel': 'VecStim', 'numCells': 10, 'spkTimes': spkt_Ad01}  # input from Adelta
netParams.popParams['C_PEP01'] = {'cellModel': 'VecStim', 'numCells': 40, 'spkTimes': spkt_C01}  # input from peptidergic C fibers
netParams.popParams['C_NP01'] = {'cellModel': 'VecStim', 'numCells': 40, 'spkTimes': spkt_C01}  # input from non-peptidergic C fibers

## SPINAL NEURONS
netParams.popParams['PKC01' ] = {'cellType': 'PKC', 'numCells': 30} # PKCg+ neurons (excitatory)
netParams.popParams['VGLUT301'] = {'cellType': 'EXdl', 'numCells': 5} # VGLUT3+ neurons (excitatory)
netParams.popParams['PV01'] = {'cellType': 'IN', 'numCells': 15} # PV+ neurons (inhibitory)
netParams.popParams['DOR01' ] = {'cellType': 'EXdl', 'numCells': 30} # DOR+ neurons (excitatory)
netParams.popParams['TrC01'] = {'cellType': 'EXib', 'numCells': 10} # Transient Central neurons (excitatory)
netParams.popParams['DYN01'] = {'cellType': 'IN', 'numCells': 60} # Central/DYN+ neurons (inhibitory)
netParams.popParams['SOM01'] = {'cellType': 'SOM', 'numCells': 15} # SOM+ neurons (excitatory)
netParams.popParams['CR01'] = {'cellType': 'CR', 'numCells': 20} # CR+ neurons (excitatory)
netParams.popParams['ISLET01'] = {'cellType': 'IN', 'numCells': 15} # Islet-type neurons (inhibitory)
netParams.popParams['NK101'] = {'cellType': 'PRO', 'numCells': 10} # NK1+ neurons (projection)

#--** NETWORK 02 **--#
## PRIMARY AFFERENTS 
netParams.popParams['Ab_SAI02'] = {'cellModel': 'VecStim', 'numCells': 5, 'spkTimes': spkt_SAI02}  # input from Ab_slow adapting type I
netParams.popParams['Ab_SAII02'] = {'cellModel': 'VecStim', 'numCells': 5, 'spkTimes': spkt_SAII02}  # input from Ab_slow adapting type I
netParams.popParams['Ad02'] = {'cellModel': 'VecStim', 'numCells': 10, 'spkTimes': spkt_Ad02}  # input from Adelta
netParams.popParams['C_PEP02'] = {'cellModel': 'VecStim', 'numCells': 40, 'spkTimes': spkt_C02}  # input from peptidergic C fibers
netParams.popParams['C_NP02'] = {'cellModel': 'VecStim', 'numCells': 40, 'spkTimes': spkt_C02}  # input from non-peptidergic C fibers

# SPINAL NEURONS
netParams.popParams['PKC02' ] = {'cellType': 'PKC02', 'numCells': 30} # PKCg+ neurons (excitatory)
netParams.popParams['VGLUT302'] = {'cellType': 'EXdl', 'numCells': 5} # VGLUT3+ neurons (excitatory)
netParams.popParams['PV02'] = {'cellType': 'IN', 'numCells': 15} # PV+ neurons (inhibitory)
netParams.popParams['DOR02' ] = {'cellType': 'EXdl', 'numCells': 30} # DOR+ neurons (excitatory)
netParams.popParams['TrC02'] = {'cellType': 'EXib', 'numCells': 10} # Transient Central neurons (excitatory)
netParams.popParams['DYN02'] = {'cellType': 'IN', 'numCells': 60} # Central/DYN+ neurons (inhibitory)
netParams.popParams['SOM02'] = {'cellType': 'SOM', 'numCells': 15} # SOM+ neurons (excitatory)
netParams.popParams['CR02'] = {'cellType': 'CR', 'numCells': 20} # CR+ neurons (excitatory)
netParams.popParams['ISLET02'] = {'cellType': 'IN', 'numCells': 15} # Islet-type neurons (inhibitory)
netParams.popParams['NK102'] = {'cellType': 'PRO', 'numCells': 10} # NK1+ neurons (projection)


###################################################################################################################################
#   Synaptic Mechanisms
###################################################################################################################################
netParams.defaultThreshold = -30

netParams.synMechParams['AMPA'] = {'mod': 'AMPA_DynSyn'   , 'tau_rise': 0.1, 'tau_decay': 5            }
netParams.synMechParams['NMDA'] = {'mod': 'NMDA_DynSyn'   , 'tau_rise': 2  , 'tau_decay': 100          }
netParams.synMechParams['NK13'] = {'mod': 'NK1_DynSyn'    , 'tau_rise': 100, 'tau_decay': 1000         }
netParams.synMechParams['NK23'] = {'mod': 'NK1_DynSyn'    , 'tau_rise': 200, 'tau_decay': 3000         }
netParams.synMechParams['GABA'] = {'mod': 'GABAa_DynSyn'  , 'tau_rise': 0.1, 'tau_decay': 20, 'e': -70 }
netParams.synMechParams['GLY']  = {'mod': 'Glycine_DynSyn', 'tau_rise': 0.1, 'tau_decay': 10, 'e': -70 }

###################################################################################################################################
#   Connect the spinal cord stimulation/natural input
###################################################################################################################################
# FOR 50mN
# from Ab_SAI to interneurons
netParams.connParams['Ab_SAI_AMPA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'PKC01'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'PKC01'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->VGLUT301'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'VGLUT301'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->VGLUT301'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'VGLUT301'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->PV01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'PV01'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->PV01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'PV01'},  
    'weight': cfg.Ab_IN_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->DYN01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI01'}, 
    'postConds': {'popLabel': 'DYN01'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

# from Ab_SAII to interneurons
netParams.connParams['Ab_SAII_AMPA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'PKC01'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'PKC01'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->VGLUT301'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'VGLUT301'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->VGLUT301'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'VGLUT301'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->PV01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'PV01'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->PV01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'PV01'},  
    'weight': cfg.Ab_IN_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->DYN01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII01'}, 
    'postConds': {'popLabel': 'DYN01'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

# from Ad to interneurons
netParams.connParams['Ad_AMPA->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'DOR01'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'DOR01'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

netParams.connParams['Ad_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

netParams.connParams['Ad_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

# from peptidergic C to interneurons
netParams.connParams['C_PEP_AMPA->TrC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'TrC01'},  
    'weight': cfg.C_TrC_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_AMPA->DYN01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'DYN01'},  
    'weight': cfg.C_DYN_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->ISLET01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'ISLET01'},  
    'weight': cfg.C_ISLET_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.C_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.C_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_NK1->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.C_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NK13'}

# from non-peptidergic C to interneurons
netParams.connParams['C_NP_AMPA->TrC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'TrC01'},  
    'weight': cfg.C_TrC_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'SOM01'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_NP_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'CR01'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_NP_AMPA->DYN01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'DYN01'},  
    'weight': cfg.C_DYN_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_AMPA->ISLET01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP01'}, 
    'postConds': {'popLabel': 'ISLET01'},  
    'weight': cfg.C_ISLET_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

###################################################################################################################################
#   Connect spinal neurons (inh>ex, ex>ex, ex/inh>projection)
###################################################################################################################################
# TO PKC NEURONS
netParams.connParams['VGLUT3_AMPA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'PKC01'},  
    'weight': cfg.VGLUT3_PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'PKC01'},  
    'weight': cfg.VGLUT3_PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PV_GABA->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV01'}, 
    'postConds': {'popLabel':'PKC01'},  
    'weight': cfg.PV_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['PV_GLY->PKC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV01'}, 
    'postConds': {'popLabel':'PKC01'},  
    'weight': cfg.PV_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

# TO DOR NEURONS
netParams.connParams['PV_GABA->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV01'}, 
    'postConds': {'popLabel':'DOR01'},  
    'weight': cfg.PV_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['PV_GLY->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV01'}, 
    'postConds': {'popLabel':'DOR01'},  
    'weight': cfg.PV_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['VGLUT3_AMPA->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'DOR01'},  
    'weight': cfg.VGLUT3_PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->DOR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'DOR01'},  
    'weight': cfg.VGLUT3_PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO TrC NEURONS
netParams.connParams['PKC_AMPA->TrC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC01'}, 
    'postConds': {'popLabel':'TrC01'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['ISLET_GABA->TrC01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'ISLET01'}, 
    'postConds': {'popLabel':'TrC01'},  
    'weight': cfg.ISLET_GABA,       
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO ISLET NEURONS
netParams.connParams['DYN_GABA->ISLET01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel':'ISLET01'},  
    'weight': cfg.DYN_ISLET_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO DYN NEURONS
netParams.connParams['ISLET_GABA->DYN01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'ISLET01'}, 
    'postConds': {'popLabel':'DYN01'},  
    'weight': cfg.ISLET_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO SOM NEURONS
netParams.connParams['DYN_GABA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.DYN_EX_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.DYN_EX_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['TrC_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.TrC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['TrC_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.TrC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PKC_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['PKC_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['VGLUT3_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.VGLUT3_SOM_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT301'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.VGLUT3_SOM_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['DOR_AMPA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.DOR_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['DOR_NMDA->SOM01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR01'}, 
    'postConds': {'popLabel':'SOM01'},  
    'weight': cfg.DOR_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO CR NEURONS
netParams.connParams['DYN_GABA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.DYN_EX_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.DYN_EX_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['TrC_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.TrC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['TrC_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.TrC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PKC_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['PKC_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['DOR_AMPA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.DOR_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['DOR_NMDA->CR01'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR01'}, 
    'postConds': {'popLabel':'CR01'},  
    'weight': cfg.DOR_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO NK1 NEURONS
netParams.connParams['SOM_AMPA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['SOM_NMDA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['SOM_NK1->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NK13'}

netParams.connParams['CR_AMPA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['CR_NMDA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['CR_NK1->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.EX_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NK13'}

netParams.connParams['DYN_GABA->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.DYN_NK1_GABA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->NK101'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN01'}, 
    'postConds': {'popLabel': 'NK101'},  
    'weight': cfg.DYN_NK1_GLY,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

# --****************--#
# --** NETWORK 02 **--#
# --****************--#
# ##################################################################################################################################
#   Connect the spinal cord stimulation/natural input
# ##################################################################################################################################
# from Ab_SAI to interneurons
netParams.connParams['Ab_SAI_AMPA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'PKC02'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'PKC02'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->VGLUT302'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'VGLUT302'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->VGLUT302'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'VGLUT302'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->PV02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'PV02'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAI_NMDA->PV02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'PV02'},  
    'weight': cfg.Ab_IN_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAI_AMPA->DYN02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAI02'}, 
    'postConds': {'popLabel': 'DYN02'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

# from Ab_SAII to interneurons
netParams.connParams['Ab_SAII_AMPA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'PKC02'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'PKC02'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->VGLUT302'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'VGLUT302'},  
    'weight': cfg.Ab_EX_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->VGLUT302'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'VGLUT302'},  
    'weight': cfg.Ab_EX_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->PV02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'PV02'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

netParams.connParams['Ab_SAII_NMDA->PV02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'PV02'},  
    'weight': cfg.Ab_IN_NMDA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'NMDA'} 

netParams.connParams['Ab_SAII_AMPA->DYN02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ab_SAII02'}, 
    'postConds': {'popLabel': 'DYN02'},  
    'weight': cfg.Ab_IN_AMPA,           
    'sec': 'dend',
    'probability': 0.2,
    'delay': 1.0,
    'loc': 0.5,
    'synMech': 'AMPA'} 

# from Ad to interneurons
netParams.connParams['Ad_AMPA->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'DOR02'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'DOR02'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

netParams.connParams['Ad_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

netParams.connParams['Ad_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.Ad_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'AMPA'}

netParams.connParams['Ad_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'Ad02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.Ad_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 5.0,
    'loc': 0.50,
    'synMech': 'NMDA'}

# from peptidergic C to interneurons
netParams.connParams['C_PEP_AMPA->TrC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'TrC02'},  
    'weight': cfg.C_TrC_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_AMPA->DYN02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'DYN02'},  
    'weight': cfg.C_DYN_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->ISLET02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'ISLET02'},  
    'weight': cfg.C_ISLET_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_AMPA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.C_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_PEP_NMDA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.C_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_PEP_NK1->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_PEP02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.C_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NK13'}

# from non-peptidergic C to interneurons
netParams.connParams['C_NP_AMPA->TrC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'TrC02'},  
    'weight': cfg.C_TrC_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'SOM02'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_NP_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.C_EX_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'CR02'},  
    'weight': cfg.C_EX_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['C_NP_AMPA->DYN02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'DYN02'},  
    'weight': cfg.C_DYN_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['C_NP_AMPA->ISLET02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'C_NP02'}, 
    'postConds': {'popLabel': 'ISLET02'},  
    'weight': cfg.C_ISLET_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 10.0,
    'loc': 0.5,
    'synMech': 'AMPA'}

##################################################################################################################################
#  Connect spinal neurons (inh>ex, ex>ex, ex/inh>projection)
##################################################################################################################################
# TO PKC02 NEURONS
netParams.connParams['VGLUT3_AMPA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'PKC02'},  
    'weight': cfg.VGLUT3_PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'PKC02'},  
    'weight': cfg.VGLUT3_PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PV_GABA->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV02'}, 
    'postConds': {'popLabel':'PKC02'},  
    'weight': cfg.PV_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['PV_GLY->PKC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV02'}, 
    'postConds': {'popLabel':'PKC02'},  
    'weight': cfg.PV_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

# TO DOR NEURONS
netParams.connParams['PV_GABA->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV02'}, 
    'postConds': {'popLabel':'DOR02'},  
    'weight': cfg.PV_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['PV_GLY->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PV02'}, 
    'postConds': {'popLabel':'DOR02'},  
    'weight': cfg.PV_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['VGLUT3_AMPA->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'DOR02'},  
    'weight': cfg.VGLUT3_PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->DOR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'DOR02'},  
    'weight': cfg.VGLUT3_PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO TrC NEURONS
netParams.connParams['PKC_AMPA->TrC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC02'}, 
    'postConds': {'popLabel':'TrC02'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['ISLET_GABA->TrC02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'ISLET02'}, 
    'postConds': {'popLabel':'TrC02'},  
    'weight': cfg.ISLET_GABA,       
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO ISLET NEURONS
netParams.connParams['DYN_GABA->ISLET02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel':'ISLET02'},  
    'weight': cfg.DYN_ISLET_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO DYN NEURONS
netParams.connParams['ISLET_GABA->DYN02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'ISLET02'}, 
    'postConds': {'popLabel':'DYN02'},  
    'weight': cfg.ISLET_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

# TO SOM NEURONS
netParams.connParams['DYN_GABA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.DYN_EX_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.DYN_EX_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['TrC_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.TrC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['TrC_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.TrC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PKC_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['PKC_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['VGLUT3_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.VGLUT3_SOM_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['VGLUT3_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'VGLUT302'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.VGLUT3_SOM_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['DOR_AMPA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.DOR_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['DOR_NMDA->SOM02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR02'}, 
    'postConds': {'popLabel':'SOM02'},  
    'weight': cfg.DOR_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO CR NEURONS
netParams.connParams['DYN_GABA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.DYN_EX_GABA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.DYN_EX_GLY,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}

netParams.connParams['TrC_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.TrC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['TrC_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'TrC02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.TrC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['PKC_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.PKC_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['PKC_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'PKC02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.PKC_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['DOR_AMPA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.DOR_AMPA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['DOR_NMDA->CR02'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DOR02'}, 
    'postConds': {'popLabel':'CR02'},  
    'weight': cfg.DOR_NMDA,           
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

# TO NK1 NEURONS
netParams.connParams['SOM_AMPA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['SOM_NMDA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['SOM_NK1->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'SOM02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NK13'}

netParams.connParams['CR_AMPA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_AMPA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'AMPA'}

netParams.connParams['CR_NMDA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_NMDA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NMDA'}

netParams.connParams['CR_NK1->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'CR02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.EX_NK1_NK1,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'NK13'}

netParams.connParams['DYN_GABA->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.DYN_NK1_GABA,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GABA'}

netParams.connParams['DYN_GLY->NK102'] = {
    'oneSynPerNetcon': True,
    'preConds': {'popLabel': 'DYN02'}, 
    'postConds': {'popLabel': 'NK102'},  
    'weight': cfg.DYN_NK1_GLY,
    'probability': 0.2,
    'sec': 'dend',
    'delay': 0.5, 
    'loc': 0.5,
    'synMech': 'GLY'}