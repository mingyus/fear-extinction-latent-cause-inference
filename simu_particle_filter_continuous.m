function simu_particle_filter_continuous(i_expList,update,pars,maxpost,N_particles,N_simu)

rng(0);
expname = {'sp','re'};
if maxpost; maxpost_str = 'maxpost'; else; maxpost_str = 'nomaxpost'; end

for i_exp = i_expList
    %% Free parameters
    if strcmp(update, 'RL')  % Recorla-Wagner learning
        
        [alpha,A,slope,baserate,eta0t,eta1t,eta0s,eta1s,v0t,v0s,rep] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7),pars(8),pars(9),pars(10),pars(11));
        eta0 = [eta0t,eta0s]; eta1 = [eta1t,eta1s]; % learning rates
        v0 = [v0t,v0s]; % the probability of tone and shock for a new cause
        
        filename = [maxpost_str '_' update '_Nparticles' num2str(N_particles) '_Nsimu' num2str(N_simu)...
            '_alpha' num2str(alpha) '_A' num2str(A) 'slope' num2str(slope) 'baserate' num2str(baserate)...
            'eta0t' num2str(eta0t) 'eta1t' num2str(eta1t) 'eta0s' num2str(eta0s) 'eta1s' num2str(eta1s)...
            'v0t' num2str(v0t) 'v0s' num2str(v0s)];
        
        % binary vs continuous shock
        if length(pars) > 11
            sigma = pars(12);
            filename = [filename 'sigma' num2str(sigma) '_' expname{i_exp}];
        else
            sigma = nan;
            filename = [filename '_' expname{i_exp}];
        end
        
    elseif strcmp(update, 'inference')  % Bayesian inference of bernoulli probability with dirichlet prior (beta)
        
        [alpha,A,slope,baserate,beta0t,beta1t,beta0s,beta1s,rep] = deal(pars(1),pars(2),pars(3),pars(4),pars(5),pars(6),pars(7),pars(8),pars(9));
        beta = [beta0t,beta1t; beta0s,beta1s];
        
        filename = [maxpost_str '_' update '_Nparticles' num2str(N_particles) '_Nsimu' num2str(N_simu)...
            '_alpha' num2str(alpha) '_A' num2str(A) 'slope' num2str(slope) 'baserate' num2str(baserate)...
            'beta0t' num2str(beta0t) 'beta1t' num2str(beta1t) 'beta0s' num2str(beta0s) 'eta1s' num2str(beta1s)...
            '_' expname{i_exp}];
    end
    
    %% Trial information: number of trials, trial indices, and timing information
    
    N_trials_train = 3;
    N_trials_extinction = 24;
    
    switch i_exp
        case 1 % spontaneous recovery
            
            % Procedure: conditioning, 24h, extinction, 24h, long term memory test, 30d, spontaneous recovery test
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
            interval_trials = [(t_tone+ITI)*ones(1, N_trials_train - 1), 24*3600, (t_tone+ITI)*ones(1,N_trials_extinction-1), ...
                24*3600, (t_tone+ITI)*ones(1,N_trials_lmtest-1), 30*24*3600, (t_tone+ITI)*ones(1,N_trials_srtest-1)]/3600;
            t_trials = [0, cumsum(interval_trials)];
            
        case 2 % reinstatement
            
            % Procedure: conditioning, 24h, extinction, 24h, reinstatement (unsignaled shocks), 24h, test
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
            interval_trials = [(t_tone+ITI)*ones(1, N_trials_train - 1), 24*3600, (t_tone+ITI)*ones(1,N_trials_extinction-1), ...
                24*3600, (t_tone+ITI)*ones(1,N_trials_reinstatement-1), 24*3600, (t_tone+ITI)*ones(1,N_trials_lmtest-1)]/3600;
            t_trials = [0, cumsum(interval_trials)];
            
    end
    
    % standard extinction, deconditioning, reduced intensity
    shocks = {
        zeros(1,N_trials_extinction),
        [-(0:N_trials_extinction/2-1)/(N_trials_extinction/2) + 1, zeros(1,N_trials_extinction/2)],
        [exp(-(0:N_trials_extinction/2-1) / (N_trials_extinction/5)), zeros(1,N_trials_extinction/2)]
    };
    shocks{4} = [shocks{2}(1:N_trials_extinction/2) + shocks{2}(N_trials_extinction/2:-1:1) - shocks{3}(N_trials_extinction/2:-1:1), zeros(1,N_trials_extinction/2)];
    
    ind_shock = {[1,3,6,10,15],[1,6,10,13,15]};
    for i_cond = 1:length(ind_shock)
        shocks{i_cond+4} = zeros(1,N_trials_extinction);
        shocks{i_cond+4}(ind_shock{i_cond}) = 1;
    end
    
    N_exp_condition = length(shocks);
    
    %% Simulation
    
    % Features: tone(f_{t,1}) and shock(f_{t,2}); both take values 0 (absent) and 1 (present)
    N_features = 2;
    N_featurevalues = 2;
    
    predict_shock_all = zeros(N_simu, N_trials, N_exp_condition);
    N_ce = N_trials_train + N_trials_extinction;
    if i_exp == 1
        N_re = 0;
    else
        N_re = N_trials_reinstatement;
    end
%     p_one_cause_pre = zeros(N_exp_condition, N_ce);
%     p_two_causes1_pre = zeros(N_exp_condition, N_ce);
%     p_two_causes2_pre = zeros(N_exp_condition, N_ce);
    p_one_cause_post = zeros(N_exp_condition, N_ce);
    p_two_causes_post = zeros(N_exp_condition, N_ce);
%     p_two_causes2_post = zeros(N_exp_condition, N_ce);
    cause_assignment = zeros(N_exp_condition, N_ce+N_re);
%     p_cause = zeros(N_exp_condition, 3);
%     p_cause_lmtest = zeros(N_exp_condition, 3);
%     p_cause_reins = zeros(N_exp_condition);
    
    for i_exp_condition = 1:N_exp_condition
        
        % SET FEATURE VALUES
        F = zeros(N_trials,2);
        % set training trials: w/ tone and shock
        F(ind_trials_train,1) = 1; F(ind_trials_train,2) = 1;
        % set extinction trials: all w/ tone, some w/ shock
        F(ind_trials_extinction,1) = 1;
        F(ind_trials_extinction,2) = shocks{i_exp_condition};
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
        
        % # particles
        particle_set = zeros(1,N_particles);
        N_samples = N_particles;
        
        for i_simu = 1:N_simu
            
            % initialize weights, causes, observations, values/counts
            c = zeros(N_trials,N_particles);
            N_causes = zeros(1,N_particles);
            predict = nan(1,N_particles); % shock prediction (per particle)
            predict_shock = zeros(1,N_trials); % shock prediction
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
                    [c(i_trial,i_particle),isnew] = rand_ddCRP(alpha,A,slope,baserate,c(1:i_trial-1,i_particle),N_causes(i_particle),t_trials(1:i_trial));
                    
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
                        
                        if ~isnan(sigma) && i_feature == 2  % only for continuous shocks
                            likelihood(i_feature,i_particle) = normpdf(F(i_trial,i_feature), pfeature, sigma);
                        else
                            likelihood(i_feature,i_particle) = pfeature * (F(i_trial,i_feature)==1) + (1-pfeature) * (F(i_trial,i_feature)==0);
                        end
                        
                        if i_feature == 2  % record prediction for shock
                            predict(i_particle) = pfeature;
                        end
                        
                        logw(i_particle) = logw(i_particle) + log(likelihood(i_feature,i_particle));                       
                        
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
                predict_shock(i_trial) = sum(likelihood(1,:)./sum(likelihood(1,:),2).*predict);
                predict_shock_baseline(i_trial) = mean(predict);
                
%                 % record cause sequence probabilities for one-cause and two-cause (after cue)
%                 p_one_cause_pre(i_exp_condition, i_trial) = mean(sum(c(1:i_trial,:) == 1, 1) == i_trial);
%                 if i_trial <= N_trials_train
%                     p_two_causes1_pre(i_exp_condition, i_trial) = 0;
%                     p_two_causes2_pre(i_exp_condition, i_trial) = 0;
%                 else
%                     idx_cause1 = ind_trials_train;
%                     idx_cause2 = (N_trials_train+1) : i_trial;
%                     p_two_causes1_pre(i_exp_condition, i_trial) = mean(sum(c(idx_cause1,:) == 1,1)==length(idx_cause1) & sum(c(idx_cause2,:) == 2,1)==length(idx_cause2));
%                     
%                     idx_cause1 = intersect([ind_trials_train ind_shock{i_exp_condition}+N_trials_train], 1:i_trial);
%                     idx_cause2 = setdiff(1:i_trial, idx_cause1);
%                     p_two_causes2_pre(i_exp_condition, i_trial) = mean(sum(c(idx_cause1,:) == 1,1)==length(idx_cause1) & sum(c(idx_cause2,:) == 2,1)==length(idx_cause2));
%                 end
%                 if i_trial == N_trials - 3 - 4 && i_exp == 1
%                     for i_cause = 1:3
%                         p_cause_lmtest(i_exp_condition, i_cause) = mean(c(i_trial,:) == i_cause);
%                     end
%                 end
                if i_trial == N_trials - 3  % the first test trial
                    cause_assignment(i_exp_condition, 1:i_trial-1) = c(1:i_trial-1,1);  % cause assignment for all trials before test
%                     for i_cause = 1:3
%                         p_cause(i_exp_condition, i_cause) = mean(c(i_trial,:) == i_cause);
%                     end
%                     if i_exp == 2
%                         p_cause_reins(i_exp_condition) = mean(c(i_trial,:) == c(i_trial-1,:));
%                     end
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
                
                % record cause sequence probabilities for one-cause and two-cause (end of trial)
                p_one_cause_post(i_exp_condition, i_trial) = mean(sum(c(1:i_trial,:) == 1, 1) == i_trial);
                if i_trial <= N_trials_train
                    p_two_causes_post(i_exp_condition, i_trial) = 0;
                else
                    % conditioning: c1, extinction: c2
                    idx_cause1 = ind_trials_train;
                    idx_cause2 = (N_trials_train+1) : i_trial;
                    p_two_causes_post(i_exp_condition, i_trial) = mean(sum(c(idx_cause1,:) == 1,1)==length(idx_cause1) & sum(c(idx_cause2,:) == 2,1)==length(idx_cause2));
                end
                
                % take the maximum of posterior (end of each session/day)
                if i_trial < N_trials && maxpost && interval_trials(i_trial) > 23
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
                
            end
            
            if i_exp == 2
                predict_shock(N_trials_train+N_trials_extinction+1:N_trials_train+N_trials_extinction+N_trials_reinstatement) = nan;
            end
            predict_shock_all(i_simu,:,i_exp_condition) = predict_shock;
            predict_shock_all_baseline(i_simu,:,i_exp_condition) = predict_shock_baseline;
            
        end
        
        % plot the overall freeze-rate prediction curves
        p_shock = mean(predict_shock_all(:,:,i_exp_condition),1);
        p_freeze = func_pshock2freeze(p_shock);
        if rep > 0
            for i_trial = 2:N_trials
                if interval_trials(i_trial - 1) < 23
                    p_freeze(i_trial) = rep * p_freeze(i_trial-1) + (1-rep) * func_pshock2freeze(p_shock(i_trial));
                end
            end
        end
        
    end
    
    %% save simulation results
    save(['results_continuous/' filename '.mat'], 'predict_shock_all', 'predict_shock_all_baseline', 'p_one_cause_post', 'p_two_causes_post', 'cause_assignment');
    %, 'p_one_cause_pre', 'p_two_causes1_pre', 'p_two_causes2_pre', 'p_one_cause_post', 'p_two_causes1_post', 'p_two_causes2_post', 'cause_assignment', 'p_cause', 'p_cause_lmtest', 'p_cause_reins'
    
end
end


function [c,isnew,pmf] = rand_ddCRP(alpha,A,slope,baserate,c_old,N_causes,t)

if isempty(c_old)
    c = 1; isnew = 1; pmf = 1;
    return;
end

% current time point
t_current = t(end);
% previous time points
t_old = t(1:end-1);

% probability for each cause
pmf = zeros(1,N_causes+1);
for i_cause = 1:N_causes
    deltat = t_current - t_old(c_old==i_cause);
    pmf(i_cause) = A*sum(exp(-slope*deltat)) + baserate;
end
pmf(end) = alpha; % new cause
pmf = pmf/sum(pmf); % normalization

% sample
cmf = cumsum(pmf);
c = find(rand()<=cmf,1,'first');

if c > N_causes
    isnew = 1;
else
    isnew = 0;
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
