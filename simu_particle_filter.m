function simu_particle_filter(i_expList,update,pars,maxpost,N_particles,N_simu,ifPlot)

rng('shuffle');
colors = [0,0,255; 61,121,4; 217,0,0]/255;
expname = {'sp','re'};

for i_exp = i_expList
    %% Free parameters
    if strcmp(update, 'RL')  % Recorla-Wagner learning
        
        [alpha,slope,baserate,eta0t,eta1t,eta0s,eta1s,v0t,v0s] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7),pars(8),pars(9));
        eta0 = [eta0t,eta0s]; eta1 = [eta1t,eta1s]; % learning rates
        v0 = [v0t,v0s]; % the probability of tone and shock for a new cause
        
        filename = [update '_Nparticles' num2str(N_particles) '_Nsimu' num2str(N_simu)...
            '_alpha' num2str(alpha) 'slope' num2str(slope) 'baserate' num2str(baserate)...
            'eta0t' num2str(eta0t) 'eta1t' num2str(eta1t) 'eta0s' num2str(eta0s) 'eta1s' num2str(eta1s)...
            'v0t' num2str(v0t) 'v0s' num2str(v0s)...
            '_' expname{i_exp}];
        
    elseif strcmp(update, 'inference')  % Bayesian inference of bernoulli probability with dirichlet prior (beta)
        
        [alpha,slope,baserate,beta0t,beta1t,beta0s,beta1s] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7));
        beta = [beta0t,beta1t; beta0s,beta1s];
        
        filename = [update '_Nparticles' num2str(N_particles) '_Nsimu' num2str(N_simu)...
            '_alpha' num2str(alpha) 'slope' num2str(slope) 'baserate' num2str(baserate)...
            'beta0t' num2str(beta0t) 'beta1t' num2str(beta1t) 'beta0s' num2str(beta0s) 'eta1s' num2str(beta1s)...
            '_' expname{i_exp}];
    end
    
    %% Trial information: number of trials, trial indices, and timing information
    
    N_exp_condition = 3;
    ind_shock = {[],[1,3,6,10,15],[1,6,10,13,15]};  % standard extinction, gradual extinction, gradual reverse
    
    switch i_exp
        case 1 % spontaneous recovery
            
            % Procedure: conditioning (3 trials), 24h, extinction (24 trials), 24h, long term memory test (4 trials), 30d, spontaneous recovery test (4 trials)
            N_trials_train = 3;
            N_trials_extinction = 24;
            N_trials_lmtest = 4;
            N_trials_srtest = 4;
            N_trials = N_trials_train + N_trials_extinction + N_trials_lmtest + N_trials_srtest;
            ind_trials_train = 1:N_trials_train;
            ind_trials_extinction = (ind_trials_train(end)+1):(ind_trials_train(end)+N_trials_extinction);
            ind_trials_lmtest = (ind_trials_extinction(end)+1):(ind_trials_extinction(end)+N_trials_lmtest);
            ind_trials_srtest = (ind_trials_lmtest(end)+1):(ind_trials_lmtest(end)+N_trials_srtest);
            session_name = cell(1,N_trials);
            session_name(ind_trials_train) = {'train'};
            session_name(ind_trials_extinction) = {'extinction'};
            session_name(ind_trials_lmtest) = {'long-term memory test'};
            session_name(ind_trials_srtest) = {'spontaneous recovery test'};
            
            t_tone = 20;
            ITI = 160;
            interval_trials = [t_tone+ITI, t_tone+(ITI+40), 24*3600, (t_tone+ITI)*ones(1,N_trials_extinction-1), ...
                24*3600, (t_tone+ITI)*ones(1,N_trials_lmtest-1), 30*24*3600, (t_tone+ITI)*ones(1,N_trials_srtest-1)]/3600;
            t_trials = [0, cumsum(interval_trials)];
            
        case 2 % reinstatement
            
            % Procedure: conditioning (3 trials), 24h, extinction (24 trials), 24h, 2 unsignaled shocks, 24h, test (4 trials)
            N_trials_train = 3;
            N_trials_extinction = 24;
            N_trials_reinstatement = 2;
            N_trials_lmtest = 4;
            N_trials = N_trials_train + N_trials_extinction + N_trials_reinstatement + N_trials_lmtest;
            ind_trials_train = 1:N_trials_train;
            ind_trials_extinction = (ind_trials_train(end)+1):(ind_trials_train(end)+N_trials_extinction);
            ind_trials_reinstatement = (ind_trials_extinction(end)+1):(ind_trials_extinction(end)+N_trials_reinstatement);
            ind_trials_lmtest = (ind_trials_reinstatement(end)+1):(ind_trials_reinstatement(end)+N_trials_lmtest);
            session_name(ind_trials_train) = {'train'};
            session_name(ind_trials_extinction) = {'extinction'};
            session_name(ind_trials_reinstatement) = {'reinstatement'};
            session_name(ind_trials_lmtest) = {'test'};
            
            t_tone = 20;
            ITI = 160;
            interval_trials = [t_tone+ITI, t_tone+(ITI+40), 24*3600, (t_tone+ITI)*ones(1,N_trials_extinction-1), ...
                24*3600, (t_tone+ITI)*ones(1,N_trials_reinstatement-1), 24*3600, (t_tone+ITI)*ones(1,N_trials_lmtest-1)]/3600;
            t_trials = [0, cumsum(interval_trials)];
            
    end
    
    %% Simulation
    
    h = figure();
    if ifPlot
        h1 = figure('pos',[100 200 1200 400]);
        h2 = figure();
    end
    
    % Features: tone(f_{t,1}) and shock(f_{t,2}); both take values 0 (absent) and 1 (present)
    N_features = 2;
    N_featurevalues = 2;
    
    for i_exp_condition = 1:N_exp_condition
        
        % SET FEATURE VALUES
        F = zeros(N_trials,2);
        % set training trials: w/ tone and shock
        F(ind_trials_train,1) = 1; F(ind_trials_train,2) = 1;
        % set extinction trials: all w/ tone, some w/ shock
        F(ind_trials_extinction,1) = 1;
        if ~isempty(ind_shock{i_exp_condition})
            F(ind_trials_extinction(ind_shock{i_exp_condition}),2) = 1;
        end
        % set reinstatement trials: all w/o tone, all w/ shock (only exp2)
        if i_exp == 2
            F(ind_trials_reinstatement,1) = 0; F(ind_trials_reinstatement,2) = 1;
        end
        % set test trials: all w/ tone, all w/o shock
        F(ind_trials_lmtest,1) = 1; F(ind_trials_lmtest,2) = 0;
        % set spontaneous recovery trials: all w/ tone, all w/o shock
        if i_exp == 1
            F(ind_trials_srtest,1) = 1; F(ind_trials_srtest,2) = 0;
        end
        % Feature values used in prediction (all trials are set to be w/ shock)
        F_predict = F;
        F_predict(:,2) = 1;
        
        % # particles
        particle_set = zeros(1,N_particles);
        N_samples = N_particles;
        
        predict_shock_all = zeros(N_simu,N_trials);
        
        for i_simu = 1:N_simu
            
            % initialize weights, causes, observations, values/counts
            c = zeros(N_trials,N_particles);
            N_causes = zeros(1,N_particles);
            predict = nan(N_features,N_particles);
            predict_shock = zeros(1,N_trials);
            logw = zeros(1,N_particles);
            if strcmp(update, 'RL')
                value = nan(N_trials,N_features,10,N_particles);
            elseif strcmp(update, 'inference')
                counts = nan(N_features,N_featurevalues,10,N_particles);
                counts_total = nan(N_features,10,N_particles);
            end
            
            for i_trial = 1:N_trials
                
                likelihood = zeros(N_features,N_particles);
                
                for i_particle = 1:N_particles
                    
                    % SAMPLE CAUSE: For each particle, sampling a hypothetical cluster assignment for the next observation from CRP
                    [c(i_trial,i_particle),isnew] = rand_ddCRP(alpha,slope,baserate,c(1:i_trial-1,i_particle),N_causes(i_particle),t_trials(1:i_trial));
                    
                    if isnew
                        N_causes(i_particle) = N_causes(i_particle)+1;
                        % initiating values/counts for this new cause
                        if strcmp(update, 'RL')
                            value(i_trial,:,c(i_trial,i_particle),i_particle) = v0;
                        elseif strcmp(update, 'inference')
                            counts(:,:,c(i_trial,i_particle),i_particle) = beta;
                            counts_total(:,c(i_trial,i_particle),i_particle) = sum(beta,2);
                            value(i_trial,:,c(i_trial,i_particle),i_particle) = beta(:,2)./sum(beta,2);
                        end
                    end
                    
                    % LIKELIHOOD: Recomputing the weights given the next observation based on the likelihood (remember to normalize)
                    for i_feature = 1:N_features
                        
                        % use values/counts to calculate likelihood for both tone and shock separately
                        if strcmp(update, 'RL')
                            pfeature = value(i_trial,i_feature,c(i_trial,i_particle),i_particle);
                        elseif strcmp(update, 'inference')
                            i_featurevalue = 1 * (F(i_trial,i_feature) == 0) + 2 * (F(i_trial,i_feature) == 1);
                            counts(i_feature,i_featurevalue,c(i_trial,i_particle),i_particle) = counts(i_feature,i_featurevalue,c(i_trial,i_particle),i_particle) + 1;
                            counts_total(i_feature,c(i_trial,i_particle),i_particle) = counts_total(i_feature,c(i_trial,i_particle),i_particle) + 1;
                            pfeature = counts(i_feature,2,c(i_trial,i_particle),i_particle) / counts_total(i_feature,c(i_trial,i_particle),i_particle);
                        end
                        
                        likelihood(i_feature,i_particle) = pfeature * (F(i_trial,i_feature)==1) + (1-pfeature) * (F(i_trial,i_feature)==0);
                        logw(i_particle) = logw(i_particle) + log(likelihood(i_feature,i_particle));
                        
                        predict(i_feature,i_particle) = pfeature * (F_predict(i_trial,i_feature)==1) + (1-pfeature) * (F_predict(i_trial,i_feature)==0);
                        
                        % RL update
                        if strcmp(update, 'RL')
                            % active cause: update values
                            eta = (F(i_trial,i_feature)==0)*eta0(i_feature) + (F(i_trial,i_feature)==1)*eta1(i_feature);
                            value(i_trial+1,i_feature,c(i_trial,i_particle),i_particle) = value(i_trial,i_feature,c(i_trial,i_particle),i_particle)...
                                + eta*(F(i_trial,i_feature) - value(i_trial,i_feature,c(i_trial,i_particle),i_particle));
                            % non-active causes: values don't change
                            all_causes = 1:max(c(1:i_trial-1,i_particle));
                            value(i_trial+1,i_feature,all_causes~=c(i_trial,i_particle),i_particle) = value(i_trial,i_feature,all_causes~=c(i_trial,i_particle),i_particle);
                        end
                    end
                end
                
                % model predicetion per trial
                predict_shock(i_trial) = sum(likelihood(1,:)./sum(likelihood(1,:),2).*predict(2,:));
                
                if ifPlot
                    figure(h1);
                    % visualization of cause probabilities (prior)
                    subplot(1,3,1);
                    plot_pcause(c(i_trial,:));
                    title('prior');
                    % visualization of tone/shock probability associated with each cause, and the overall likelihood
                    subplot(1,3,2);
                    plot_likelihood(likelihood,c(i_trial,:),exp(logw));
                    title('likelihood');
                end
                
                % RESAMPLING:
                % normalizing w
                if sum(exp(logsumexp(logw))) == 0
                    logw = - log(N_particles) * ones(N_particles);
                else
                    logw = logw - logsumexp(logw);
                end
                
                % Sampling (with replacement) from the current set of particles according to the importance weights
                if i_trial>1
                    rnds = rand(1,N_samples);
                    for i_sample = 1:N_samples
                        particle_set(i_sample) = find(rnds(i_sample)<=cumsum(exp(logw)),1,'first');
                    end
                    c = c(:,particle_set);
                    N_causes = N_causes(particle_set);
                    if strcmp(update, 'RL')
                        value = value(:,:,:,particle_set);
                    elseif strcmp(update, 'inference')
                        counts = counts(:,:,:,particle_set);
                        counts_total = counts_total(:,:,particle_set);
                    end
                end
                
                logw = zeros(1,N_particles); % re-initializing w after resampling
                
                % take the maximum of posterior (end of each session/day)
                if (maxpost && interval_trials(i_trial) > 23)
                    % find the cause sequence with the highest posterior
                    allcauseseq = unique(c(1:i_trial,:)','rows');
                    numcause = zeros(1,size(allcauseseq,1));
                    idx_cause = zeros(1,N_particles);
                    exampleparticle = zeros(1,size(allcauseseq,1));
                    for icause = 1:size(allcauseseq,1)
                        for i_particle = 1:N_particles
                            if sum(allcauseseq(icause,:) ~= c(1:i_trial,i_particle)')==0
                                numcause(icause) = numcause(icause) + 1;
                                idx_cause(i_particle) = icause;
                                this_particle = i_particle;
                            end
                        end
                        exampleparticle(icause) = this_particle;
                    end
                    N_maxpost = length(find(numcause == max(numcause)));
                    particle_maxpost = exampleparticle(numcause == max(numcause));
                    p_sample = 1 / N_maxpost * ones(1,N_maxpost);
                    rnds = rand(1,N_samples);
                    for i_sample = 1:N_samples
                        particle_set(i_sample) = particle_maxpost(find(rnds(i_sample)<=cumsum(p_sample),1,'first'));
                    end
                    c = c(:, particle_set);
                    N_causes = N_causes(particle_set);
                    if strcmp(update, 'RL')
                        value = value(:,:,:,particle_set);
                    elseif strcmp(update, 'inference')
                        counts = counts(:,:,:,particle_set);
                        counts_total = counts_total(:,:,particle_set);
                    end
                end
                
                
                if ifPlot
                    % visualization of cause probabilities (posterior)
                    subplot(1,3,3);
                    plot_pcause(c(i_trial,:));
                    title('posterior');
                    sgtitle([session_name{i_trial} ': ' num2str(F(i_trial,2))],'FontSize',15);
                    % visulization of posterior of cause sequence
                    figure(h2);
                    plot_pcauseseq(c(1:i_trial,:));
                    pause();
                end
                
            end
            
            if i_exp == 2
                predict_shock(28:29) = nan;
            end
            predict_shock_all(i_simu,:) = predict_shock;
            
        end
        
        % plot the overall freeze-rate prediction curves
        figure(h); hold on;
        p(i_exp_condition) = plot(1:N_trials,mean(predict_shock_all,1),'-o','linewidth',1.5,'color',colors(i_exp_condition,:));
        
    end
    
    %% set format for freeze-rate prediction plot
    figure(h);
    ylim([0 1])
    xlabel('Trial index')
    ylabel('Estimated probability of shock')
    set(gca,'fontsize',15)
    colors = [0,0,255; 61,121,4; 217,0,0]/255;
    if i_exp == 1
        scatter([1:N_trials_train],0.05*ones(1,N_trials_train),[],colors(1,:),'filled');
        scatter([1:N_trials_train N_trials_train+ind_shock{2}],0.1*ones(1,N_trials_train+length(ind_shock{2})),[],colors(2,:),'filled');
        scatter([1:N_trials_train N_trials_train+ind_shock{3}],0.15*ones(1,N_trials_train+length(ind_shock{3})),[],colors(3,:),'filled');
    else
        ind_unsignaledShock = N_trials_train+N_trials_extinction+1:N_trials_train+N_trials_extinction+2;
        scatter([1:N_trials_train ind_unsignaledShock],0.05*ones(1,N_trials_train+2),[],colors(1,:),'filled');
        scatter([1:N_trials_train N_trials_train+ind_shock{2} ind_unsignaledShock],0.1*ones(1,N_trials_train+length(ind_shock{2})+2),[],colors(2,:),'filled');
        scatter([1:N_trials_train N_trials_train+ind_shock{3} ind_unsignaledShock],0.15*ones(1,N_trials_train+length(ind_shock{3})+2),[],colors(3,:),'filled');
    end
    
    changeday = find(interval_trials>=24)+1;
    for i = 1:numel(changeday)
        line([changeday(i)-0.5 changeday(i)-0.5],[0,1],'color',[0.5 0.5 0.5]);
    end
    
    legend(p,{'standard extinction','gradual extinction','gradual reverse'})
    legend boxoff;
    
    print(['figure/' filename '.png'],'-dpng');
end
end

function s = logsumexp(X, dim)
% Compute log(sum(exp(X),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Mo Chen (sth4nth@gmail.com).
if nargin == 1
    % Determine which dimension sum will use
    dim = find(size(X)~=1,1);
    if isempty(dim), dim = 1; end
end
% subtract the largest in each dim
y = max(X,[],dim);
s = y+log(sum(exp(bsxfun(@minus,X,y)),dim));
i = isinf(y);
if any(i(:))
    s(i) = y(i);
end
end

function plot_pcause(c)
[n_c,~] = histcounts(c,0.5:1:max(c+0.5));
bar(1:max(c),n_c/sum(n_c),'FaceColor',[0.5 0.5 0.5],'BarWidth',0.5);
xlim([0.5 max(c+0.5)]); ylim([0 1]);
set(gca,'ytick',0:0.1:1);
xlabel('Cause'); ylabel('Probability');
set(gca,'FontSize',15);
end

function plot_likelihood(likelihood,c,w)
p_tone = zeros(1,max(c)); p_shock = zeros(1,max(c)); p_both = zeros(1,max(c));
for i_cause = 1:max(c)
    if isempty((c==i_cause))
        p_tone(i_cause) = nan;
        p_shock(i_cause) = nan;
        p_both(i_cause) = nan;
    else
        p_tone(i_cause) = mean(likelihood(1,c==i_cause));
        p_shock(i_cause) = mean(likelihood(2,c==i_cause));
        p_both(i_cause) = mean(w(c==i_cause));
    end
end
if max(c) == 1
    bar(1:max(c)+1,[p_tone,p_shock,p_both;nan,nan,nan]); set(gca,'xtick',1);
else
    bar(1:max(c),[p_tone;p_shock;p_both]');
end
xlim([0.5 max(c+0.5)]); ylim([0 1]);
xlabel('Cause'); ylabel('Probability');
set(gca,'FontSize',15);
end

function plot_pcauseseq(c)
clf; hold on;
N_particles = size(c,2);
allcauseseq = unique(c','rows');
numcause = zeros(1,size(allcauseseq,1));
legend_text = cell(1,size(allcauseseq,1));
i = 0;
for icause = 1:size(allcauseseq,1)
    for i_particle = 1:N_particles
        if sum(allcauseseq(icause,:) ~= c(:,i_particle)')==0
            numcause(icause) = numcause(icause) + 1;
        end
    end
    if numcause(icause) > N_particles/100
        i = i+1;
        bar(i,numcause(icause),'FaceColor',[0.5 0.5 0.5],'BarWidth',0.5);
        legend_text{i} = [num2str(i) ': ' num2str(allcauseseq(icause,:))];
    end
end
if i > 0
    xlim([0.5 i+0.5]); ylim([0 N_particles]);
    set(gca,'ytick',0:N_particles/10:N_particles,'yticklabel',0:0.1:1,'xtick',1:i);
    legend(legend_text(1:i));
end
xlabel('Cause sequence'); ylabel('Probability');
set(gca,'FontSize',15);
end