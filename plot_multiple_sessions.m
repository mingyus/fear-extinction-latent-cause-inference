clear; close all;

simuname = 'multiple_sessions';
expName = {'sp', 're'};
exp_type = {'Spontaneous recovery', 'Reinstatement'};

Nsubplots = 2;
h = figure('Position', [0,0,600,600]);
fontsize1 = 12;
fontsize2 = 18;

% exp info
Nc = 3;
Ne = 9;
Nr = 2;
Nt = 4;
N_trials = {[3,9,4], [3,9,2,4]};

marker = {'x', '^'};

for isubplot = 1:Nsubplots
    %% load data
    rep = 0.7;  % perseveration probability
    
    iExp = isubplot;
    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.05_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', simuname,'_', expName{iExp}, '.mat']);
    
    
    if iExp == 1 % remove memory test trials
        predict_shock_all = predict_shock_all(:,[1:Nc+Ne, Nc+Ne+Nt+1:end],:);
    end
    
    p_shock = mean(predict_shock_all,1);
    p_freeze = func_pshock2freeze(p_shock);
    if rep > 0
        for i_trial = 2:size(p_freeze,2)
            if (iExp == 1 && ~ismember(i_trial, [Nc+1,Nc+Ne+1])) || (iExp == 2 && ~ismember(i_trial, [Nc+1,Nc+Ne+1,Nc+Ne+Nr+1]))
                p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
            end
        end
    end
    
    
    %% figure
    figure(h); subplot(Nsubplots,1,isubplot); hold on; colors = get(gca,'colororder');
    
    for iCond = 1:size(p_freeze,3)
        start = 1; x_trials_all = [];
        for day = 1:length(N_trials{iExp})
            finish = start + N_trials{iExp}(day) - 1;
            trial_idx = start:finish;
            x_trials = (start:finish) + day;
            x_trials_all = [x_trials_all, x_trials];
            if iCond == 2 && day == 2
                for subday = 1:3
                    p(iCond) = plot(x_trials((subday-1)*3+1:subday*3), mean(p_freeze(:,trial_idx((subday-1)*3+1:subday*3),iCond),1),'-o','linewidth',1.5,'color',colors(2+iCond,:));
                end
                
            else
                p(iCond) = plot(x_trials, mean(p_freeze(:,trial_idx,iCond),1),'-o','linewidth',1.5,'color',colors(2+iCond,:));
            end
            start = finish + 1;
        end
    end
    
    changeday = cumsum(N_trials{iExp}(1:end-1)) + (1:length(N_trials{iExp})-1) + 1.5;
    if iExp == 1
        text([1.2, Nc+Ne/2+1.2, Nc+Ne+Nt/2+3], -0.2*ones(1,3), {'Conditioning', 'Extinction', 'Test'}, 'fontsize', fontsize2);
        xticks(x_trials_all)
        xticklabels({'1','2','3','1','2','3','4','5','6','7','8','9','1','2','3','4'})
    else
        text([1, Nc+Ne/2+1.2, Nc+Ne+Nr/2+2, Nc+Ne+Nr+Nt/2+4], [-0.15, -0.15, -0.192, -0.15]-0.05, {'Conditioning', 'Extinction', compose('Reinstate-\n   ment'), 'Test'}, 'fontsize', fontsize2);
        xticks(x_trials_all)
        xticklabels({'1','2','3','1','2','3','4','5','6','7','8','9','1','2','1','2','3','4'})
    end
    
    for i = 1:numel(changeday)
        line([changeday(i)-0.5 changeday(i)-0.5],[0,0.9],'color',[0.5 0.5 0.5]*(~ismember(i,[2,3])) + [0.75 0.75 0.75]*ismember(i,[2,3]),'linestyle','--');
    end
    
    if isubplot == 1
        legend(p, {'Massed extinction','Spaced extinction'},'Position',[0.72 0.9 0.1 0.08])
        legend boxoff;
    end
    
%     title(exp_type{iExp});
    
    ylim([0 1]);
    xlim([0, sum(N_trials{iExp})+length(N_trials{iExp})+1]);
    ylabel('Predicted freezing rate')
    
    yh = get(gca,'ylabel');
    pos = get(yh,'position');
    pos(1) = pos(1) - 1;
    set(yh,'position',pos);
    
    set(gca,'fontsize', fontsize2);
    
end