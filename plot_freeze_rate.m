function plot_freeze_rate(model)
expName = {'sp', 're'};
exp_str = {'Spontaneous recovery', 'Reinstatement'};
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
fontsize1 = 12;
fontsize2 = 18;

for isubplot = 1:Nsubplots
    %% load data
    rep = 0.7;  % perseveration probability
    
    switch model
        case 'main'
            iExp = isubplot;
            load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
            exp_type = {'Spontaneous recovery', 'Reinstatement'};
        case 'CNN'
            iExp = isubplot;
            load('results/CNN_p_freeze_tone.mat');
            exp_type = {'Spontaneous recovery', 'Reinstatement'};
        case 'alternative'
            iExp = 1; % plot only SP
            iCond = isubplot;
            switch iCond
                case 1
                    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'Main model';
                case 2
                    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'Standard CRP prior';
                case 3
                    load(['results/maxpost_inference_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1beta0t0.5beta1t0.5beta0s0.95eta1s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'Learning through inference';
                case 4
                    load(['results/nomaxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'Full posterior';
                case 5
                    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'No perseveration';
                    rep = 0;
                case 6
                    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.1eta1t0.1eta0s0.1eta1s0.2v0t0.5v0s0.05_', expName{iExp}, '.mat']);
                    model_variant = 'No perseveration + low learning rate';
                    rep = 0;
            end
    end
    
    if ismember(model, {'main', 'alternative'})
        p_shock = mean(predict_shock_all,1);
        p_freeze = func_pshock2freeze(p_shock);
        p_shock_baseline = mean(predict_shock_all_baseline,1);
        p_freeze_baseline = func_pshock2freeze(p_shock_baseline);
        if rep > 0
            for i_trial = 2:size(p_freeze,2)
                if (iExp == 1 && ~ismember(i_trial, [4,28,32])) || (iExp == 2 && ~ismember(i_trial, [4,30]))
                    p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
                    p_freeze_baseline(:,i_trial, :) = rep * p_freeze_baseline(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
                end
            end
        end
    end
    
    %% figure
    switch model
        case 'main'
            figure(h); subplot(Nsubplots,1,isubplot); hold on;
        case 'alternative'
            figure(h(isubplot)); plt = subplot(1,1,1); hold on;
            plt.Position = plt.Position .* [1 1 1 0.8] + [0 0.1 0 0];
        case 'CNN'
            figure(h); subplot(1,Nsubplots,isubplot); hold on;
    end
    
    if ismember(model, {'main', 'alternative'})
        
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
            if strcmp(model, 'main') || (strcmp(model, 'alternative') && ismember(isubplot, [3,6]))
                text([-0.8, 14, 29.5, 35.5], -0.15*ones(1,4), {'Conditioning', 'Extinction', 'Memory', 'Test'}, 'fontsize', fontsize2);
            end
            xticks(x_trials{iExp})
            xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','3','4','1','2','3','4'})
        else
            changeday = [4+0.5, 28+1.5, 30+2.5];
            if strcmp(model, 'main') || (strcmp(model, 'alternative') && ismember(isubplot, [3,6]))
                text([-0.8, 13, 28, 33.5], [-0.15, -0.15, -0.192, -0.15], {'Conditioning', 'Extinction', compose('Reinstate-\n   ment'), 'Test'}, 'fontsize', fontsize2);
            end
            xticks(x_trials{iExp})
            xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','1','2','3','4'})
        end
        
        for i = 1:numel(changeday)
            line([changeday(i)-0.5 changeday(i)-0.5],[0,0.9],'color',[0.5 0.5 0.5],'linestyle','--');
        end
        
        if isubplot == 1
            switch model
                case 'main'
                    legend(p, {'Standard extinction','Gradual extinction','Gradual reverse'},'Position',[0.8 0.9 0.1 0.08])
                case 'alternative'
                    legend(p, {'Standard extinction','Gradual extinction','Gradual reverse'},'Position',[0.76 0.85 0.1 0.08])
            end
            legend boxoff;
        end
        
        switch model
            case 'main'
                title(exp_type{iExp});
            case 'alternative'
                title(model_variant);
        end
        
        ylim([0 1]);
        xlim([0, N_trials_total(iExp)+4]);
        ylabel('Predicted freezing rate')
        
        yh = get(gca,'ylabel');
        pos = get(yh,'position');
        pos(1) = pos(1) - 1.5;
        set(yh,'position',pos);
        
        set(gca,'fontsize', fontsize2);
        
    else
        ind_shock = {[],[1,3,6,10,15],[1,6,10,13,15]};
        N_trials_extinction = size(p_freeze_all{1, 1},1);
        
        for iCond = 1:3
            x = (1:N_trials_extinction)';
            y = mean(p_freeze_all{iExp, iCond},2);
            err = std(p_freeze_all{iExp, iCond},0,2)/sqrt(NRats);
            f(iCond) = plot(x, y, '-o', 'linewidth', 1.5, 'color', colors(iCond,:));
            fill([x;flip(x)], [y-err;flip(y+err)], colors(iCond,:), 'LineStyle','none', 'facealpha', 0.2);
        end
        
        scatter(ind_shock{2}, 0.1*ones(1,5), [], colors(2,:),'filled');
        scatter(ind_shock{3}, 0.15*ones(1,5), [], colors(3,:),'filled');
        
        if iExp == 1
            legend(f, {'Standard extinction','Gradual extinction','Gradual reverse'},'Position',[0.35 0.73 0.1 0.08])
            legend boxoff
        end
        
        ylim([0 1]);
        xlim([0, N_trials_extinction+1]);
        ylabel('Proportion of time freezing')
        xlabel('Extinction trial index');
        title(exp_type{iExp});
        
        yh = get(gca,'ylabel');
        pos = get(yh,'position');
        pos(1) = pos(1) - 1.3;
        set(yh,'position',pos);
        
        set(gca,'fontsize', fontsize2);
        
    end
    
end
end