from pyjmi import transfer_optimization_problem
import matplotlib.pyplot as plt
import numpy as N
from pymodelica import compile_fmu
from pyfmi import load_fmu


reference = compile_fmu("StartupOptimizationAtt.StartupAttReference",
                 ("CombinedCycle.mo","StartupOptimizationAtt.mop"),
                 target='me',compiler='optimica',
                 compiler_options={'enable_variable_scaling':False})

model_fmu = load_fmu(reference)
res_fmu = model_fmu.simulate(start_time=0.,final_time=10000.)


#plot simulation
plant_p = res_fmu['plant.p']
plant_sigma = res_fmu['plant.sigma']
plant_load = res_fmu['plant.load']
plant_TITs = res_fmu['plant.TIT']
plant_TITm = res_fmu['plant.turbine.TIT']
plant_w_att = res_fmu['plant.w_att']
time = res_fmu['time']

plt.figure(1)
plt.clf()
plt.subplot(3,1,1)
plt.plot(time,plant_p/1e6)
plt.ylabel('p [MPa]')
plt.grid(True)
plt.subplot(3,1,2)
plt.plot(time,plant_sigma/1e6)
plt.grid(True)
plt.ylabel('stress [MPa]')
plt.subplot(3,1,3)
plt.plot(time,plant_load)
plt.grid(True)
plt.ylabel('load [normalized]')
plt.xlabel('time [s]')

plt.figure(2)
plt.clf()
plt.plot(time,plant_TITs,time,plant_TITm)
plt.legend(('set point','measured'))
plt.grid(True)
plt.ylabel('Set point and measured values of vapour temperature at turbine inlet [K]')
plt.xlabel('time [s]')

plt.figure(3)
plt.clf()
plt.plot(time,plant_w_att)
plt.grid(True)
plt.ylabel('Attemperation Mass Flow [kg/s]')
plt.xlabel('time [s]')



#Optimization
compiler_options={'equation_sorting':True}
optim = transfer_optimization_problem("StartupOptimizationAtt.StartupAtt",("CombinedCycle.mo","StartupOptimizationAtt.mop"), compiler_options)
#option setting
n_e = 40
n_cp = 1
opt_options = optim.optimize_options()
opt_options['init_traj'] = res_fmu
opt_options['n_e'] = n_e
opt_options['n_cp'] = n_cp
opt_options['blocking_factors'] = N.ones(n_e, dtype=int)
#opt_options['result_mode'] = 'element_interpolation'
#opt_options['IPOPT_options']['mu_strategy'] = 'monotone'
#opt_options['IPOPT_options']['mumps_mem_percent'] = 10000

# variable elimination
eliminables = optim.getEliminableVariables()
algebraics = optim.getVariables(optim.REAL_ALGEBRAIC)
print "Number of algebraics:  ",len(algebraics)
eliminable_algebraics = [a for a in algebraics if a in eliminables]
variables_to_eliminate = list()
variables_with_bounds = list()
for v in eliminable_algebraics:
    if not v.hasAttributeSet("min") and not v.hasAttributeSet("max"):
        variables_to_eliminate.append(v)
    else:
        variables_with_bounds.append(v)
# Elimination of unbounded variables
optim.markVariablesForElimination(variables_to_eliminate)
# Elimination of bounded variables (skip elimination of plant.sigma)
bounded_eliminations = [v for v in variables_with_bounds if not v.getName()=="plant.sigma"]
optim.markVariablesForElimination(bounded_eliminations)
optim.eliminateVariables()
eliminated_vars = optim.getEliminatedVariables()
print "Number of variables that were eliminated: ", len(eliminated_vars)

#solution
res_opt = optim.optimize(options=opt_options)


#plot simulation
plant_p = res_opt['plant.p']
plant_sigma = res_opt['plant.sigma']
plant_load = res_opt['plant.load']
plant_TITs = res_opt['plant.TIT']
plant_TITm = res_opt['plant.turbine.TIT']
plant_w_att = res_opt['plant.w_att']
turbine_Tvap = res_opt['plant.turbine.wat_vap.T']
turbine_Tint = res_opt['plant.turbineShaft.T[1]']
turbine_Tm = res_opt['plant.turbineShaft.Tm']
time = res_opt['time']


plt.figure(4)
plt.clf()
plt.subplot(3,1,1)
plt.plot(time,plant_p/1e6)
plt.ylabel('p [MPa]')
plt.grid(True)
plt.subplot(3,1,2)
plt.plot(time,plant_sigma/1e6)
plt.legend(('Attemperation','Comparison'))
plt.grid(True)
plt.ylabel('stress [MPa]')
plt.subplot(3,1,3)
plt.plot(time,plant_load)
plt.grid(True)
plt.ylabel('load [normalized]')
plt.xlabel('time [s]')

plt.figure(5)
plt.clf()
plt.plot(time,plant_TITs,time,plant_TITm)
plt.legend(('set point','measured'))
plt.grid(True)
plt.ylabel('Set point and measured values of vapour temperature at turbine inlet [K]')
plt.xlabel('time [s]')

plt.figure(6)
plt.clf()
plt.plot(time,plant_w_att)
plt.grid(True)
plt.ylabel('Attemperation Mass Flow [kg/s]')
plt.xlabel('time [s]')

plt.figure(7)
plt.clf()
plt.plot(time,turbine_Tvap, time, turbine_Tint, time, turbine_Tm)
plt.legend(('vapour','internal boundary', 'mean'))
plt.ylabel('Rotor Temperatures [K]')
plt.xlabel('time [s]')
plt.grid(True)
plt.show()