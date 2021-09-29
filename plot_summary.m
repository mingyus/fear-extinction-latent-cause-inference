clear; close all;

expName = {'sp', 're'};
exp_str = {'Spontaneous recovery', 'Reinstatement'};
cond_str = {'Standard extinction','Gradual extinction','Gradual reverse'};
colors = [0,0,255; 61,121,4; 217,0,0]/255;

h1 = figure('Position', [200,200,800,500]);
h2 = figure('Position', [200,200,800,500]);

fontsize1 = 12;
fontsize2 = 18;

for iExp = 1:2
    %% load data
    % experiment
    switch iExp
        case 1
            data = xlsread('handscores/gershman2013_experiment1.xlsx');
        case 2
            data = xlsread('handscores/gershman2013_experiment2.xlsx');
    end
    
    % model simulation
    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_', expName{iExp}, '.mat']);
    rep = 0.7;  % perseveration probability
    p_shock = mean(predict_shock_all,1);
    p_freeze = func_pshock2freeze(p_shock);
    p_shock_baseline = mean(predict_shock_all_baseline,1);
    p_freeze_baseline = func_pshock2freeze(p_shock_baseline);
    
    for i_trial = 2:size(p_freeze,2)
        if (iExp == 1 && ~ismember(i_trial, [4,28,32])) || (iExp == 2 && ~ismember(i_trial, [4,30]))
            p_freeze(:,i_trial, :) = rep * p_freeze(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
            p_freeze_baseline(:,i_trial, :) = rep * p_freeze_baseline(:,i_trial-1, :) + (1-rep) * func_pshock2freeze(p_shock(:,i_trial, :));
        end
    end
    
    for iType = 1:2
        for iCond = 1:3
            if iType == 2  % replication of original paper
                if iExp == 1
                    idxTrial = {2:5, 7:10, 21, 22:25};
                    idxRat = {1:16, 22:37, [43:49 51:58]};
                else
                    idxTrial = {2:5, 7:10, 13, 14:17};
                    idxRat = {35:42, 1:12, 18:29};
                end
            else  % model simulation
                if iExp == 1
                    idxTrial = {4:7, 24:27, 32, 32:35};
                else
                    idxTrial = {4:7, 24:27, 30, 30:33};
                end
            end
        end

        
        %% summary of p(shock) curve
        figure(h1); subplot(2,2,(iExp-1)*2+iType); hold on;
        for iCond = 1:3
            x = {1:4, 6:9, 11, 12:15};
            if iType == 2  % replication of original paper
                for iSession = 1:4
                    y = nanmean(data(idxRat{iCond},idxTrial{iSession}),1);
                    sem = nanstd(data(idxRat{iCond},idxTrial{iSession}),0,1)./sqrt(sum(~isnan(data(idxRat{iCond},idxTrial{iSession})),1));
                    if iSession == 3
                        f(iCond) = errorbar(x{iSession}, y, sem, '-o', 'linewidth', 1.5, 'color', colors(iCond,:));
                    else
                        f(iCond) = plot(x{iSession}, y, '-o', 'linewidth', 1.5, 'color', colors(iCond,:));
                        fill([x{iSession},flip(x{iSession})], [y-sem,flip(y+sem)], colors(iCond,:), 'LineStyle','none', 'facealpha', 0.2);
                    end
                end
            else  % model simulation
                for iSession = 1:4
                    if iSession == 3
                        y = mean(p_freeze_baseline(:, idxTrial{iSession}, iCond),1);
                    else
                        y = mean(p_freeze(:, idxTrial{iSession}, iCond),1);
                    end
                    f(iCond) = plot(x{iSession}, y, '-o', 'linewidth', 1.5, 'color', colors(iCond,:));
                end
            end
        end
        
        plot([5, 5], [0,1], '--', 'color', [0.5,0.5,0.5]);
        plot([10, 10], [0,1], '--', 'color', [0.5,0.5,0.5]);
        xlim([0,16]); ylim([0,1]);
        xticks(1:15); xticklabels({'1','2','3','4','','21','22','23','24','','0','1','2','3','4'})
        yticks(0:0.2:1); yticklabels(0:20:100)
        if iExp == 1 && iType == 2
            legend(f, cond_str,'Position',[0.83 0.83 0.1 0.08]);
            legend boxoff
        end
        ylabel('Freezing (%)');
        xlabel('     Extinction                         Test');
        set(gca,'fontsize',fontsize1);
        title(exp_str{iExp});
        
        %% average of four trials of test - last four trials of extinction
        figure(h2); subplot(2,2,(iExp-1)*2+iType); hold on;
        if iType == 2  % replication of original paper
            y = []; sem = [];
            for iCond = 1:3
                diff = nanmean(data(idxRat{iCond}, idxTrial{4}) - data(idxRat{iCond}, idxTrial{2}),2);
                y(iCond) = nanmean(diff, 1);
                sem(iCond) = nanstd(diff, 0, 1)./sqrt(sum(~(isnan(data(idxRat{iCond},idxTrial{2}(1)))|isnan(data(idxRat{iCond},idxTrial{4}(1)))),1));
            end
            errorbar(1:3, y, sem, 'k', 'linestyle', 'none', 'linewidth', 1.5);
        else  % model simulation
            y = squeeze(mean(mean(p_freeze(:, idxTrial{4}, :),1) - mean(p_freeze(:, idxTrial{2}, :),1), 2));
        end
        bar(1:3, y, 'FaceColor', 'k', 'FaceAlpha', 0.5);
        xticks([1,2,3])
        my_xticklabels(gca, 1:3, {compose('Standard         \nextinction         '),compose('Gradual\nextinction'),compose('         Gradual\n         reverse')}, 'fontsize', fontsize1);
        xlim([0.2,3.8]);
        ylim([-0.22,0.52]);
        ylabel('Test - Ext: Freezing (%)');
        title(exp_str{iExp});
        set(gca,'fontsize',fontsize1);
        
    end
end