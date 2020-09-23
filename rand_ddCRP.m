function [c,isnew,pmf] = rand_ddCRP(slope,alpha,baserate,c_old,N_causes,t)

if isempty(c_old)
    c = 1; isnew = 1; pmf = 1;
    return;
end

% current time point
t_current = t(end);
% previous time points
t_old = t(1:end-1);
% total number of causes so far: N_causes

% probability for each cause
pmf = zeros(1,N_causes+1);
for i_cause = 1:N_causes
    deltat = t_current - t_old(c_old==i_cause);
    pmf(i_cause) = sum(exp(-slope*deltat)) + baserate; %  + baserate*exp(-slopelong*min(deltat))
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