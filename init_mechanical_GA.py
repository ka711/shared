'''this is the code for the DHN model involved with mechanical pain network
created by K. Sekiguchi (25thMay20)
'''
from netpyne import sim
from neuron import h
					
simConfig, netParams = sim.readCmdLineArgs(simConfigDefault='cfg_mechanical_GA.py', netParamsDefault='netParams_mechanical_GA.py')

# Create network and run simulation
# sim.createSimulate(netParams = netParams, simConfig = simConfig)
sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)
