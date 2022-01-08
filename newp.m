function newp(task)

% new parameter values to test
p_original = [0.2,1,0.1,0.1,0.2,0.2,0.2,0.4,0.5,0.05,0.7];

p_new = {
    {[1], [0.1,0.15,0.2,0.3,0.5]},
    {[3], [0.05,0.07,0.1,0.2,0.5]},
    {[4], [0,0.05,0.1,0.15,0.2]},
    {[5,6,7], [0.16,0.18,0.2,0.22,0.24]},
    {[8], [0.36,0.38,0.4,0.42,0.44]},
    {[9], [0.1,0.4,0.5,0.6,0.7]},
    {[10], [0,0.02,0.05,0.08,0.1]},
    {[11], [0.5,0.6,0.7,0.8,0.9]},
    };

p = cell(1,40);
count = 0;
for i = 1:length(p_new)
    for p_change = p_new{i}{2}
        count = count + 1;
        p{count} = p_original;
        for i_p = p_new{i}{1}
            p{count}(i_p) = p_change;
        end
    end
end

if nargin > 0 && strcmp(task, 'simulation')  % simulation
    for i = 1:length(p_new)
        simu_particle_filter([1,2],'RL',p{i},1,10000,1);
    end
    
elseif nargin > 0 && strcmp(task, 'plot')  % plot
    
    par_names = {'\alpha','k','b','\eta','\eta_{shock}', 'V_0(tone)', 'V_0(shock)', 'p_r'};
    
    p0 = 0.2;
    count = 0;
    for i = 1:length(p_new)
        for p_change = p_new{i}{2}
            count = count + 1;
            plot_freeze_rate_newp('main', [p{count}, p0], [par_names{i}, ' = ', num2str(p_change)]);
            plot_cause_probability_newp(p{count}, [par_names{i}, ' = ', num2str(p_change)]);
        end
        pause(); close all;
    end
    
    par_name = 'p_0';
    for p0 = 0:0.1:0.4
        plot_freeze_rate_newp('main', [p_original, p0], [par_name, ' = ', num2str(p0)])
        plot_cause_probability_newp(p_original, [par_name, ' = ', num2str(p_change)]);
    end
    
end

end

%% plot_freeze_rate_newp
function plot_freeze_rate_newp(model, pars, p_name_value)
expName = {'sp', 're'};
colors = [0,0,255; 61,121,4; 217,0,0]/255;

switch model
    case 'main'
        Nsubplots = 2;
        h = figure('Position', [0,0,800,800]);
    case 'alternative'
        Nsubplots = 6;
        for isubplot = 1:Nsubplots
            h(isubplot) = figure('Position', [200,200,800,400]);
        end
    case 'CNN'
        Nsubplots = 2;
        h = figure('Position', [0,200,1400,400]);
end
fontsize = 18;

for isubplot = 1:Nsubplots
    % load data
    
    iExp = isubplot;
    
    [alpha,A,slope,baserate,eta0t,eta1t,eta0s,eta1s,v0t,v0s,rep,p0] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7),pars(8),pars(9),pars(10),pars(11),pars(12));
    filename = ['maxpost_RL_Nparticles10000_Nsimu1' ...
        '_alpha' num2str(alpha) '_A' num2str(A) 'slope' num2str(slope) 'baserate' num2str(baserate)...
        'eta0t' num2str(eta0t) 'eta1t' num2str(eta1t) 'eta0s' num2str(eta0s) 'eta1s' num2str(eta1s)...
        'v0t' num2str(v0t) 'v0s' num2str(v0s) '_' expName{iExp}];
    
    load(['results/' filename '.mat']);
    
    p_shock = mean(predict_shock_all,1);
    p_freeze = func_pshock2freeze(p_shock, p0);
    p_shock_baseline = mean(predict_shock_all_baseline,1);
    p_freeze_baseline = func_pshock2freeze(p_shock_baseline);
    if rep > 0
        for i_trial = 2:size(p_freeze,2)
            if (iExp == 1 && ~ismember(i_trial, [4,28,32])) || (iExp == 2 && ~ismember(i_trial, [4,30]))
                p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :),p0);
                p_freeze_baseline(:,i_trial, :) = rep * p_freeze_baseline(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :),p0);
            end
        end
    end
    
    % figure
    figure(h); subplot(Nsubplots,1,isubplot); hold on;
    
    
    N_trials_total = [35, 33];
    
    trial_idx = [{1:3, 4:27, 28:31, 32:35}; {1:3, 4:27, 28:29, 30:33}];
    x_trials = {[1:3, (4:27)+1, (28:31)+2, (32:35)+3], [1:3, (4:27)+1, (28:29)+2, (30:33)+3]};
    for iCond = 1:3
        for iphase = 1:4
            p(iCond) = plot(trial_idx{iExp,iphase}+iphase-1,mean(p_freeze(:,trial_idx{iExp,iphase},iCond),1),'-o','linewidth',1.5,'color',colors(iCond,:));
        end
    end
    
    if iExp == 1
        scatter(1:3, 0.05*ones(1,3),[],colors(1,:),'filled');
        scatter([1, 2, 3, [4, 6, 9, 13, 18]+1], 0.1*ones(1,8), [], colors(2,:),'filled');
        scatter([1, 2, 3, [4, 9, 13, 16, 18]+1], 0.15*ones(1,8), [], colors(3,:),'filled');
    else
        scatter([1, 2, 3, [28, 29]+2], 0.05*ones(1,5),[],colors(1,:),'filled');
        scatter([1, 2, 3, [4, 6, 9, 13, 18]+1, [28, 29]+2], 0.1*ones(1,10),[],colors(2,:),'filled');
        scatter([1, 2, 3, [4, 9, 13, 16, 18]+1, [28, 29]+2], 0.15*ones(1,10),[],colors(3,:),'filled');
    end
    
    if iExp == 1
        changeday = [4+0.5, 28+1.5, 32+2.5];
        if strcmp(p_name_value(1:3),'p_r')
            text([-0.8, 14, 29.5, 35.5], -0.15*ones(1,4), {'Conditioning', 'Extinction', 'Memory', 'Test'}, 'fontsize', fontsize);
        end
        xticks(x_trials{iExp})
        xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','3','4','1','2','3','4'})
    else
        changeday = [4+0.5, 28+1.5, 30+2.5];
        xticks(x_trials{iExp})
        xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','1','2','3','4'})
    end
    
    for i = 1:numel(changeday)
        line([changeday(i)-0.5 changeday(i)-0.5],[0,0.9],'color',[0.5 0.5 0.5],'linestyle','--');
    end
    
    if isubplot == 1
        title(p_name_value, 'interpreter', 'tex');
    end
    
    ylim([0 1]);
    xlim([0, N_trials_total(iExp)+4]);
    ylabel('Predicted freezing rate')
    
    yh = get(gca,'ylabel');
    pos = get(yh,'position');
    pos(1) = pos(1) - 1.5;
    set(yh,'position',pos);
    
    set(gca,'fontsize', fontsize);
    
end
end

%% plot_cause_probability_newp
function plot_cause_probability_newp(pars, p_name_value)

[alpha,A,slope,baserate,eta0t,eta1t,eta0s,eta1s,v0t,v0s,rep] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7),pars(8),pars(9),pars(10),pars(11));
filename = ['maxpost_RL_Nparticles10000_Nsimu1' ...
    '_alpha' num2str(alpha) '_A' num2str(A) 'slope' num2str(slope) 'baserate' num2str(baserate)...
    'eta0t' num2str(eta0t) 'eta1t' num2str(eta1t) 'eta0s' num2str(eta0s) 'eta1s' num2str(eta1s)...
    'v0t' num2str(v0t) 'v0s' num2str(v0s) '_sp'];
load(['results/' filename '.mat']);

colors = [0,0,255; 61,121,4; 217,0,0]/255;

p_one_cause = [p_one_cause_post(:,4:27)];
p_two_causes = [p_two_causes2_post(:,4:27)];

% end of extinction session

figure('Position', [0,0,100,150]); hold on;

for i_exp_condition = 1:3
    bar(i_exp_condition, log(p_two_causes(i_exp_condition, end) ./ p_one_cause(i_exp_condition, end)), 'FaceColor', colors(i_exp_condition,:), 'FaceAlpha', 0.85, 'edgecolor', 'none');
end

xlim([0,4]);
xticks([]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.XAxis.TickLength = [0 0];

ylim([-4,4])
ax.YAxis.Visible = 'off';
% ylabel('Log posterior ratio')

% set(gca,'fontsize',13);

set(gcf, 'color', 'none');
set(gca, 'color', 'none');

end

