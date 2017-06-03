function [ sys ] = default_opts_fade( sys )
T = size(sys.y,2);
    % Order of the Autoregressive Model
    if ~isfield(sys,'Order'); sys.Order = 2;                        end
    % Minimum number of multiplicative updates
    if ~isfield(sys,'min_iters'); sys.min_iters = 100;              end
    % Maximum number of multiplicative updates
    if ~isfield(sys,'num_iters'); sys.num_iters = 200;              end
    % Variance of noise
    if ~isfield(sys,'noise'); sys.noise = estimate_noise(sys.y);    end
    % Regularization parameter
    if ~isfield(sys,'lambda'); sys.lambda = 0;                      end
	% (p,q)-pair in the l_pq regularization of spikes
    if ~isfield(sys,'p_norm'); sys.p_norm = 1;            end
    if ~isfield(sys,'q_norm'); sys.q_norm = 1;            end
if sys.lambda > 0
    message = sprintf('Penalizing ell_%0.1f,%0.1f norm of spikes',sys.p_norm,sys.q_norm);
else
    message = 'Warning: Lambda set to 0, solving the ML problem';
end
    disp(message)
    % scale regularization parameter with respect to the length of data (from theory of lasso)
    sys.lambda = sys.lambda*sqrt(log(T)/T)*diag(sys.noise);
    % Estimate the baseline (it can also be set to be
    sys.baseline = estimate_base(sys.y,sys.noise);
    % Subtract the baseline for deconvolution
    sys.y = sys.y - sys.baseline;
    % treat everything under the baseline as noise and set it to zero
    sys.y(sys.y<0) = 0;
    sys.b0 = zeros(size(sys.y));
    % estimate the time constant nd corresponding filter parameters
    if isfield(sys,'tau')
        sys.Order = 1;
        sys.theta = [1 -exp(-1/sys.tau)];
    else
        if ~isfield(sys,'theta'); yt = sys.y'; 
            if size(sys.y,1) == 1
                sys.theta = nanmean(aryule(yt,sys.Order));
            else
            sys.theta = nanmean(aryule(yt(:,1:min(100,end)),sys.Order));   
            end
        end
    sys.tau = 1/log(-1/sum(sys.theta(2:end)));
    end
    
        sys.s0 = filter(sys.theta,1,sys.y,[],2); sys.s0(sys.s0<=0) = 1;
%     sys.s0 = ones(size(sys.y));
    
    sys.T = size(sys.y,2);
end

function sigma = estimate_noise(y)
T = size(y,2);
range_ff = [0.25; 0.5];
[pyy,ff]=pwelch(y',round(T/8),[],1000,1);
ind = ff > range_ff(1) & ff < range_ff(2);
sigma = sqrt(exp(mean(log(pyy(ind,:)/2))))';
end

function baseline = estimate_base(y,noise)
T = size(y,2);
m = repmat(min(y,[],2),1,T);
y = y-m;
noise = repmat(noise,1,T);
mask = y < 4.5*noise;
y(~mask) = nan;
baseline = nanmean(y,2);
baseline = repmat(baseline,1,T)+m;
end