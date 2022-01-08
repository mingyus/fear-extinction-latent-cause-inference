clear; close all;

filename = {'re', 'low_intensity_re'};
exp_type = {'Reinstatement', 'Reinstatement: low intensity'};
colors = [0,0,255; 61,121,4; 217,0,0]/255;

Nsubplots = 2;
h = figure('Position', [0,0,800,800]);

fontsize = 18;

for isubplot = 1:Nsubplots
    %% load data
    rep = 0.7;  % perseveration probability
    
    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', filename{isubplot}, '.mat']);
    
    p_shock = mean(predict_shock_all,1);
    p_freeze = func_pshock2freeze(p_shock);
    p_shock_baseline = mean(predict_shock_all_baseline,1);
    p_freeze_baseline = func_pshock2freeze(p_shock_baseline);
    if rep > 0
        for i_trial = 2:size(p_freeze,2)
            if ~ismember(i_trial, [4,30])
                p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
                p_freeze_baseline(:,i_trial, :) = rep * p_freeze_baseline(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
            end
        end
    end
    
    %% figure
    
    figure(h); subplot(Nsubplots,1,isubplot); hold on;
    
    N_trials_total = 33;
    
    trial_idx = {1:3, 4:27, 28:29, 30:33};
    x_trials = [1:3, (4:27)+1, (28:29)+2, (30:33)+3];
    for iCond = 1:3
        for iphase = 1:4
            p(iCond) = plot(trial_idx{iphase}+iphase-1,mean(p_freeze(:,trial_idx{iphase},iCond),1),'-o','linewidth',1.5,'color',colors(iCond,:));
        end
    end
    
    if isubplot == 1
        scatter([1, 2, 3, [28, 29]+2], 0.05*ones(1,5),[],colors(1,:),'filled');
        scatter([1, 2, 3, [4, 6, 9, 13, 18]+1, [28, 29]+2], 0.1*ones(1,10),[],colors(2,:),'filled');
        scatter([1, 2, 3, [4, 9, 13, 16, 18]+1, [28, 29]+2], 0.15*ones(1,10),[],colors(3,:),'filled');
    else
        scatter([1, 2, 3], 0.05*ones(1,3),[],colors(1,:),'filled');
        scatter([1, 2, 3, [4, 6, 9, 13, 18]+1], 0.1*ones(1,8),[],colors(2,:),'filled');
        scatter([1, 2, 3, [4, 9, 13, 16, 18]+1], 0.15*ones(1,8),[],colors(3,:),'filled');
        scatter([28, 29]+2, 0.05*ones(1,2),[],colors(1,:));
        scatter([28, 29]+2, 0.1*ones(1,2),[],colors(2,:));
        scatter([28, 29]+2, 0.15*ones(1,2),[],colors(3,:));
    end
    
    changeday = [4+0.5, 28+1.5, 30+2.5];
    text([-0.8, 13, 28, 33.5], [-0.15, -0.15, -0.192, -0.15], {'Conditioning', 'Extinction', compose('Reinstate-\n   ment'), 'Test'}, 'fontsize', fontsize);
    xticks(x_trials)
    xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','1','2','3','4'})
    
    
    for i = 1:numel(changeday)
        line([changeday(i)-0.5 changeday(i)-0.5],[0,0.9],'color',[0.5 0.5 0.5],'linestyle','--');
    end
    
    if isubplot == 2
        legend(p, {'Standard extinction','Gradual extinction','Gradual reverse'},'Position',[0.8 0.4 0.1 0.08])
        legend boxoff;
    end
    
    
    title(exp_type{isubplot});
    
    ylim([0 1]);
    xlim([0, N_trials_total+4]);
    ylabel('Predicted freezing rate')
    
    yh = get(gca,'ylabel');
    pos = get(yh,'position');
    pos(1) = pos(1) - 1.5;
    set(yh,'position',pos);
    
    set(gca,'fontsize', fontsize);
    
end