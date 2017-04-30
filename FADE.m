function sys_fade = FADE( sys )

%% An implementation of the FAst DEconvolution (FADE) Algorithm
% Input: Structure called sys with the following fields:
%        y: observed calcium traces where the first dimension is the number of ROI's and
%        second dimension is time and
%        p_norm and q_norm: Penalizes the ell_pq norm of spikes, default: p = q = 1
%        lambda: regularization parameter, default: lambda = 0 (ML Estimation)
%        Order: Autoregressive model order, defualt: Order = 2
%        theta: AR parameters, default: estimated using Yule-Walker equations
%        noise: estimate of the noise standard variation, default:
%        estimate using power spectral density
%        min_iters: minimum number of iterations, default: min_iters = 100
%        num_iters: maximum number of iterations, default: num_iters = 200
% Output: A structure containing the following fields
%        spikes: the deconvolved spikes, post processing methods should be
%        smoothed_traces: smoothed calcium traces
%        used on the spikes in most cases
%        ds: relative changes in the spikes        
%%        
sys = default_opts_fade(sys);
[p,T] = size(sys.y);

b = sys.b0;
s = sys.s0;
ds = zeros(1,sys.num_iters);
iter = 1;
ds(1) = 1;
    while iter < sys.num_iters && ( ds(iter)>0.005 || iter<sys.min_iters)
        [Lp,Ln] = calc_L(sys,s,b); % positive and negative log-likelihood terms
        Pp = calc_penalty(sys,s); % Penalty terms
        snew = Ln./(Lp+sys.lambda*Pp).*s; % multiplicative updates
        ds(iter+1) = max(abs((sum(s,2)-sum(snew,2))./sum(s,2)));
        s = snew;
        iter = iter + 1;
    end
X = filter(1,sys.theta,s,[],2);

sys_fade = sys;
sys_fade.spikes = s;
sys_fade.smoothed_traces = X + sys.baseline;
sys_fade.theta = sys.theta;
sys_fade.ds = ds(1:iter-1);
end

function [Lp,Ln] = calc_L(sys,s,b)
Ln = 2*fliplr(filter(1,sys.theta,fliplr(sys.y),[],2));
X = filter(1,sys.theta,s,[],2);
Lp = 2*fliplr(filter(1,sys.theta,fliplr(X+b),[],2));
end

function Pp = calc_penalty(sys,s)
    A1 = sum(sum(s.^sys.p_norm,2).^(sys.q_norm/sys.p_norm))^(1/sys.q_norm-1);
    A2 = repmat(sum(s.^sys.p_norm,2).^(sys.q_norm/sys.p_norm-1),1,sys.T);
    A3 = s.^(sys.p_norm-1);
    Pp = A1*A2.*A3;
end










