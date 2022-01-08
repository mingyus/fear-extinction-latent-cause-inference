clear; close all;

colors = [0,0,255; 61,121,4; 217,0,0; 0,0,0]/255;

load('results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_re.mat');
p_one_cause = p_one_cause_post(:,27);
p_two_causes = p_two_causes2_post(:,27);

load('results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_occasional_reinforcement.mat');
p_one_cause = [p_one_cause; p_one_cause_post(:,27)];
p_two_causes = [p_two_causes; p_two_causes2_post(:,27)];

figure('Position', [0,200,450,200]); hold on;

pr = log(p_two_causes(4:end) ./ p_one_cause(4:end));
violinplot(pr, 1, 'ViolinColor', [0 0 0], 'ViolinAlpha', 0.2, 'ShowData', false);

for i_exp_condition = 1:3
    bar(i_exp_condition-3, log(p_two_causes(i_exp_condition) ./ p_one_cause(i_exp_condition)), 'FaceColor', colors(i_exp_condition,:), 'FaceAlpha', 0.85, 'linestyle', 'none');
end

ylabel('Log posterior ratio')
xticks((1:4)-3)
my_xticklabels(gca, (1.03:0.95:3.88)-3.3, {compose('Standard\nextinction'),compose('Gradual\nextinction'),compose('Gradual\nreverse'),compose('Occasional\nreinforcement')}, 'fontsize', 13);

xlim([0,5]-3);
ylim([-4.5,2.5])
yticks(-4:4);
set(gca,'fontsize',13);

ax = gca;
ax.XRuler.Axle.Visible = 'off';
ax.XAxis.TickLength = [0 0];