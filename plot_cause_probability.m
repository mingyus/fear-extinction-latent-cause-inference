%% extinction: log posterior ratio (with shocks marked)

clear;
load('results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_re.mat');
colors = [0,0,255; 61,121,4; 217,0,0]/255;

p_one_cause = [p_one_cause_post(:,4:27)]; %p_one_cause_pre(:,4), 
p_two_causes = [p_two_causes2_post(:,4:27)]; %p_two_causes2_pre(:,4), 

ind_shock = {[],[1,3,6,10,15],[1,6,10,13,15]};

figure(); hold on;
for i_exp_condition = 1:3
    pl(i_exp_condition) = plot(1:24, log(p_two_causes(i_exp_condition,1:end) ./ p_one_cause(i_exp_condition,1:end)), '-o','linewidth',1.5, 'color',colors(i_exp_condition,:));
end
x = -2:30;
plot(x, zeros(length(x)), '--', 'color', 'k','linewidth',1.25)

scatter(ind_shock{2}, -2*ones(1,5), [], colors(2,:),'filled');
scatter(ind_shock{3}, -1.5*ones(1,5), [], colors(3,:),'filled');

legend(pl, {'Standard extinction','Gradual extinction','Gradual reverse'}, 'Location', 'southeast');
legend boxoff;
lpos = get(gca,'legend'); % handle to the label object
p = get(lpos,'position'); % get the current position property
p(1) = p(1) + 0.025;
p(2) = p(2) + 0.1;
set(lpos,'position',p);

ylabel('Log posterior ratio')
xlabel('Extinction trial index');

xlim([0,25]);
ylim([-2.5,2.5])
yticks(-2:2);
set(gca,'fontsize',15);

% padding of xlabel and ylabel
xh = get(gca,'xlabel');
p = get(xh,'position');
p(2) = p(2)-0.05;
set(xh,'position',p);

yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = p(1) - 0.5;
set(yh,'position',p);

%% first trial in test: cause probability
clear;

exp = {'sp', 're'};
exp_str = {'Memory', 'Test', 'Test'};

causes = {1:3, 1:2, 1:3};
colors_cause = [255, 192, 0; 165, 165, 165; 148, 32, 147]/256;

cause_labels = {
    compose('cond /\nshock'), compose('ext /\nno-shock');
    'cond + ext', 'reins';
    'shock', 'no-shock';
};

for i_plot = 1:3
    
    i_exp = 1 * (i_plot<=2) + 2 * (i_plot==3);
    
    load(['results/maxpost_RL_Nparticles10000_Nsimu1_alpha0.2_A1slope0.1baserate0.1eta0t0.2eta1t0.2eta0s0.2eta1s0.4v0t0.5v0s0.05_' exp{i_exp} '.mat']);
    
    if i_plot == 1  % long-term memory test in spontaneous recovery
       p_cause = p_cause_lmtest;
       i_exp = 1;
    end
    
    figure('Position', [200,200,300,600]);
    for i_exp_condition = 1:3
        subplot(3,1,i_exp_condition); hold on;
        barwidth = 0.7;
        if i_plot ~= 3 && i_exp_condition == 2
            barwidth = barwidth/3*2;
        end
        for i_cause = causes{i_exp_condition}
            if i_plot == 3 && i_exp_condition == 2 && i_cause == 2
                bar(i_cause, p_cause_reins(i_exp_condition), 'FaceColor', colors_cause(1,:), 'FaceAlpha', 0.3);
                bar(i_cause + 1, p_cause(i_exp_condition, i_cause), barwidth, 'FaceColor', colors_cause(i_cause+1*(i_exp_condition==2),:), 'FaceAlpha', 0.3);
                i_cause = i_cause + 1;
            else
                bar(i_cause, p_cause(i_exp_condition, i_cause), barwidth, 'FaceColor', colors_cause(i_cause+1*(i_exp_condition==2),:), 'FaceAlpha', 0.3);
            end
            
        end
        xlim([0.3, i_cause+0.7])
        xticks(1:i_cause);
        labels = {cause_labels{i_exp_condition, 1:i_cause-1}, 'new'};
        my_xticklabels(gca,1:i_cause,labels,'fontsize',15);
        ylim([0,1]);
        ylabel('Prior probability');
        if i_exp_condition == 3
            xlabel('Cause');
            xh = get(gca,'xlabel');
            p = get(xh,'position');
            p(2) = p(2)-0.15;
            set(xh,'position',p);
        end
        set(gca,'fontsize',15);
        
    end
    sgtitle(exp_str{i_plot},'fontsize',18);
end
