This repository contains data and code that accompany the manuscript entitled "Explaining the effectiveness of fear extinction through latent-cause inference".

**Table of Contents**

* `handscores/` contains the data reported in Gershman et al. (2013)
	* `FreqExtinction_Handscores_CEJ.xlsx`: Experiment 1 (spontaneous recovery)
	* `Reinstatement_FreqExt_Handscores.xlsx`: Experiment 2 (reinstatement)
* `*.m` files: code for model simulation and plotting figures (see details below)
* `results/` contains the model simulation results; can be reproduced with `simu_particle_filter.m`; used for generating figures in the manuscript.
* `stats.ipynb`: statistical test results

**Model simulation:**

`simu_particle_filter.m` contains the model simulation code. You can run the scripts below to simulate the main model and alternative models; alternatively, you can directly access the simulation results saved in `results/`.

Main model:

```
for i=1:2
	simu_particle_filter(i,'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);
end
```

Alternative models:

```
for i=1:2
	% Standard CRP prior
	simu_particle_filter(i,'RL',[0.2,1,0,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],1,10000,1);
	% Learning through inference
	simu_particle_filter(i,'inference',[0.2,1,0.1,0.1,0.5,0.5,0.95,0.05,0.7],1,10000,1);
	% Full posterior
	simu_particle_filter(i,'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7],0,10000,1);
	% No perseveration
	simu_particle_filter(i,'RL',[0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0],1,10000,1);
	% No perseveration + low learning rate
	simu_particle_filter(i,'RL',[0.2,1,0.1,0.1,0.1,0.1,0.1,0.2,0.5,0.05,0],1,10000,1,0);
end
```

**To reproduce figures in the manuscript**

* `plot_cause_probability.m`: Figure 3
* `plot_freeze_rate('main')`: Figure 4
* `plot_freeze_rate('alternative')`: Figure 6
<!-- * `plot_freeze_rate('CNN')`: Supplemenary Figure 1 --> 
* `plot_summary.m`: Figure 5, Supplemenary Figure 1

**Helper functions**

* `func_pshock2freeze.m`: the mapping from p(shock) to freezing rate.
* `my_xticklabels.m`: helper function to format xticklabels in figures (by Pekka Kumpulainen)
