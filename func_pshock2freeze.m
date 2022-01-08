function freeze = func_pshock2freeze(pshock, p0)
if nargin < 2
    p0 = 0.2;
end
% linear with offset
freeze = pshock*(1-p0)+p0;
end