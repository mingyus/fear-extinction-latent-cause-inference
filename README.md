This repository contains data and code that accompany the manuscript entitled "Explaining the effectiveness of fear extinction through latent-cause inference".

**Table of Contents**

* `handscores/` contains the data reported in Gershman et al. (2013)
	* `gershman2013_experiment1.xlsx`: Experiment 1 (spontaneous recovery)
	* `gershman2013_experiment2.xlsx`: Experiment 2 (reinstatement)
* `*.m` files: code for model simulation and plotting figures (see details below)
* `results/` contains the model simulation results; can be reproduced with `simu_particle_filter.m`; used for generating figures in the manuscript.
* `stats.ipynb`: statistical test results

**Model simulation:**

`simu_particle_filter.m` contains the model simulation code. You can run the scripts below to simulate the main model and alternative models; alternatively, you can directly access the simulation results saved in `results/`.

Main model:

```
simu_particle_filter([1,2],'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);
```

Alternative models:

```
% Standard CRP prior
simu_particle_filter([1,2],'RL',[0.2,1,0,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);

% Learning through inference
simu_particle_filter([1,2],'inference',[0.2,1,0.1,0.1,0.5,0.5,0.95,0.05,0.7],1,10000,1);

% Full posterior
simu_particle_filter([1,2],'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],0,10000,1);

% No perseveration
simu_particle_filter([1,2],'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0],1,10000,1);

% No perseveration + low learning rate
simu_particle_filter([1,2],'RL',[0.2,1,0.1,0.1,0.1,0.1,0.1,0.2,0.5,0.05,0],1,10000,1,0);
```

Other experiments:

```
% spaced vs massed extinction
simu_particle_filter_other_exps('multiple_sessions',[1,2],'RL',[0.05,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);
	
% occasional reinforcement experiments
simu_particle_filter_other_exps('occasional_reinforcement',3,'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);

% low intensity shock during reinforcement experiment
simu_particle_filter_other_exps('low_intensity',2,'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);
```

Shock with continuous intensity:

```
simu_particle_filter_continuous([1,2],'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7,0.4],1,10000,1);
```

Simulation with varying parameter values:

```
newp('simulation');
```

**To reproduce figures in the manuscript**

* `plot_cause_probability.m`: Figure 3
* `plot_freeze_rate('main')`: Figure 4
* `plot_summary.m`: Figure 5, Supplemenary Figure 2
* `plot_freeze_rate('alternative')`: Figure 6
* `plot_multiple_sessions.m`: Figure 7A
* `plot_low_intensity.m`: Figure 7B
* `plot_occasional_reinforcement.m`: Figure 7C
* `plot_continuous.m`: Figure 7D
* `newp('plot')`: Supplementary Figure 1

**Helper functions**

* `func_pshock2freeze.m`: the mapping from p(shock) to freezing rate.
* `my_xticklabels.m`: helper function to format xticklabels in figures (by Pekka Kumpulainen)
