# This is the code to create SOM positive neuron in lamina II having initial burst firing
# created by K. Sekiguchi (30th Apr 2020)

from netpyne import specs

cfg = specs.SimConfig()		# object of class SimConfig to store simulation configuration

# COMMON
cfg.Ra = 93.59353718147344           # 80
cfg.cm = 1.247406111           # 1.0
cfg.g_pas = 1.1e-05    # 1.1e-05
cfg.e_pas = -70        # -70
cfg.ek_KDRI = -84
cfg.alpha_shift_B_Na = -4.88077338 # 0
cfg.beta_shift_B_Na  = 14.65623395 # 0

## SOMA
cfg.gnabar_B_Na_soma = 0.001192007	  # 0.008
cfg.gkbar_KDRI_soma  = 0.048559067 # 0.0043
cfg.gkbar_B_DR_soma  = 0.0    # 0

## DEND
cfg.gkbar_KDRI_dend  = 0.205862692   # 0.034

## HILLOCK
cfg.gnabar_B_Na_hill = 4.954161853  # 0.73 
cfg.gkbar_B_A_hill   = 0.0   # 0
cfg.gkbar_B_DR_hill  = 0.0   # 0
cfg.gkbar_KDRI_hill  = 0.08848351852046527 # 0.076
cfg.gkbar_KDR_hill   = 0.0   # 0


cfg.hParams = {'celsius': 23, 'v_init': -81, 'alpha_shift_B_Na': cfg.alpha_shift_B_Na, 'beta_shift_B_Na': cfg.beta_shift_B_Na }

cfg.duration = 1300    # Duration of the simulation, in ms
cfg.dt = 0.025 		   # Internal integration timestep to use
cfg.verbose = False  			# Show detailed messages 
# cfg.recordTraces = {'V_soma':{'sec':'soma', 'loc': 0.5, 'var': 'v'},
#                     'gkHH2' :{'sec':'soma', 'loc': 0.5, 'var':'gk_HH2'},         #
#                     'gkBA'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_B_A'},         #
#                     'gkCa'  :{'sec':'soma', 'loc': 0.5, 'var':'gk_iKCa'},        #
#                     'gkdr'  :{'sec':'soma', 'loc': 0.5, 'var':'gkdr_borgkdr'}}   #

# cfg.recordTraces = {'V_soma':{'sec':'soma', 'loc': 0.5, 'var': 'v'},
#                     'ikHH2' :{'sec':'soma', 'loc': 0.5, 'var':'ik_HH2'},
#                     'ikBA'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_B_A'},
#                     'ikCa'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_iKCa'},
#                     'ikdr'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_borgkdr'}}

cfg.recordTraces = {'V_soma':{'sec':'soma', 'loc': 0.5, 'var': 'v'}}
                    # 'V_dend':{'sec':'dend', 'loc': 0.5, 'var': 'v'},
                    # 'V_hillock':{'sec':'hillock', 'loc': 0.5, 'var': 'v'}}
                    # 'V_axon':{'sec':'axon', 'loc': 0.5, 'var': 'v'}}
                    # 'ikHH2' :{'sec':'soma', 'loc': 0.5, 'var':'ik_HH2'},
                    # 'ikBA'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_B_A'},
                    # 'ikCa'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_iKCa'},
                    # 'ikdr'  :{'sec':'soma', 'loc': 0.5, 'var':'ik_borgkdr'}}

cfg.recordStep = 0.025 			# Step size in ms to save data (eg. V traces, LFP, etc)
cfg.filename = 'model_output'  # Set file output name
cfg.saveJson = True
cfg.savePickle = False 		# Save params, network and sim output to pickle file

cfg.analysis['plotTraces'] = {'include': ['all'], 'showFig': True, 'timeRange': [0, cfg.duration], 'oneFigPer': 'trace', 'saveFig': True}