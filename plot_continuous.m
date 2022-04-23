clear; close all;

expName = {'sp', 're'};
exp_str = {'Spontaneous recovery', 'Reinstatement'};

Nsubplots = 2;
h1 = figure('Position', [0,0,800,800]);
h2 = figure('Position', [800,400,400,300]);
fontsize = 18;

colors = [0,0,1; 0.929,0.694,0.125; 0.494,0.184,0.556; 0.6*ones(1,3); 61/255,121/255,4/255; 217/255,0,0];

linewidth = 1.5*ones(1,6);
linestyle = {'-o','-o','-o','-o','-','-'};
alpha = 0.25;

for isubplot = 1:Nsubplots
    %% exp info
    Nc = 3;
    Ne = 24;
    Nr = 2;
    Nt = 4;
    N_trials_total = [Nc+Ne+Nt, Nc+Ne+Nr+Nt];
    
    %% load data
    rep = 0.7;  % perseveration probability
    
    iExp = isubplot;
    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05sigma0.43', '_', expName{iExp}, '.mat']);
    
    if iExp == 1  % remove memory trials
        predict_shock_all = predict_shock_all(:,[1:Nc+Ne,Nc+Ne+Nt+1:end],:);
    end
    
    p_shock = mean(predict_shock_all,1);
    p_freeze = func_pshock2freeze(p_shock);
    p_shock_baseline = mean(predict_shock_all_baseline,1);
    p_freeze_baseline = func_pshock2freeze(p_shock_baseline);
    if rep > 0
        for i_trial = 2:size(p_freeze,2)
            if (iExp == 1 && ~ismember(i_trial, [Nc+1,Nc+Ne+1,Nc+Ne+Nt+1])) || (iExp == 2 && ~ismember(i_trial, [Nc+1,Nc+Ne+Nr+1]))
                p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
                p_freeze_baseline(:,i_trial, :) = rep * p_freeze_baseline(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
            end
        end
    end
    
    
    %% shock intensity
    figure(h2); hold on;
    
    shocks = {
        zeros(1,Ne),
        [-(0:Ne/2-1)/(Ne/2) + 1, zeros(1,Ne/2)],
        [exp(-(0:Ne/2-1) / (Ne/5)), zeros(1,Ne/2)]
        };
    shocks{4} = [shocks{2}(1:Ne/2) + shocks{2}(Ne/2:-1:1) - shocks{3}(Ne/2:-1:1), zeros(1,Ne/2)];
    
    for iCond = 1:size(shocks)
        plot(1:Ne,shocks{iCond},linestyle{iCond},'linewidth',linewidth(iCond),'color',colors(iCond,:));
    end
    
    ind_shock = {[1,3,6,10,15],[1,6,10,13,15]};
    for iCond = 1:length(ind_shock)
        scatter(ind_shock{iCond},ones(1,length(ind_shock{iCond}))-0.02+iCond*0.02,'o','linewidth',linewidth(iCond),'markeredgecolor',colors(iCond+4,:),'markeredgealpha',alpha);
    end
    
    yticks([0,0.5,1])
    ylim([-0.05 1.05]);
    xlim([0, Ne + 1]);
    title('Shock intensity');
    xlabel('Extinction trial index');

    set(gca,'fontsize', fontsize);
    
    %% freeze rate
    figure(h1); subplot(Nsubplots,1,isubplot); hold on;
    
    trial_idx = {{1:Nc, Nc+1:Nc+Ne, Nc+Ne+1:Nc+Ne+Nt};
        {1:Nc, Nc+1:Nc+Ne, Nc+Ne+1:Nc+Ne+Nr, Nc+Ne+Nr+1:Nc+Ne+Nr+Nt}};
    x_trials = {[1:Nc, (Nc+1:Nc+Ne)+1, (Nc+Ne+1:Nc+Ne+Nt)+2],
        [1:Nc, (Nc+1:Nc+Ne)+1, (Nc+Ne+1:Nc+Ne+Nr)+2, (Nc+Ne+Nr+1:Nc+Ne+Nr+Nt)+3]};
    for iCond = 1:size(p_freeze,3)
        for iphase = 1:length(trial_idx{iExp})
            x = trial_idx{iExp}{iphase}+iphase-1;
            y = mean(p_freeze(:,trial_idx{iExp}{iphase},iCond),1);
            p(iCond) = plot(x,y,linestyle{iCond},'linewidth',linewidth(iCond),'color',colors(iCond,:));
            if iCond > 4
                p(iCond).Color(4) = alpha;
                scatter(x,y,'o','MarkerEdgeColor',colors(iCond,:),'linewidth',linewidth(iCond),'MarkerEdgeAlpha',alpha);
            end
        end
    end
    
    if iExp == 1
        changeday = [Nc+1+0.5, Nc+Ne+1+1.5];
        text([-0.8, Nc+Ne/2-0.8, Nc+Ne+Nt/2+1.8], -0.15*ones(1,3), {'Conditioning', 'Extinction', 'Test'}, 'fontsize', fontsize);
        xticks(x_trials{iExp})
        xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','3','4'})
    else
        changeday = [Nc+1+0.5, Nc+Ne+1+1.5, Nc+Ne+Nr+1+2.5];
        text([-0.8, Nc+Ne/2-0.8, Nc+Ne+Nt/2-0.5, Nc+Ne+Nt+Nr/2+1.5], [-0.15, -0.15, -0.192, -0.15], {'Conditioning', 'Extinction', compose('Reinstate-\n   ment'), 'Test'}, 'fontsize', fontsize);
        xticks(x_trials{iExp})
        xticklabels({'1','2','3','1','','','','5','','','','','10','','','','','15','','','','','20','','','','24','1','2','1','2','3','4'})
    end
    
    for i = 1:numel(changeday)
        line([changeday(i)-0.5 changeday(i)-0.5],[0,0.9],'color',[0.5 0.5 0.5],'linestyle','--');
    end
    
%     title(exp_str{iExp});
    
    ylim([0 1]);
    xlim([0, N_trials_total(iExp)+4]);
    ylabel('Predicted freezing rate')
    
    yh = get(gca,'ylabel');
    pos = get(yh,'position');
    pos(1) = pos(1) - 1.5;
    set(yh,'position',pos);
    
    set(gca,'fontsize', fontsize);
    
end
