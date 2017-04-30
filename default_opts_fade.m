function [ sys ] = default_opts_fade( sys )
T = size(sys.y,2);
    if ~isfield(sys,'Order'); sys.Order = 2;                        end
    if ~isfield(sys,'min_iters'); sys.min_iters = 100;              end
    if ~isfield(sys,'num_iters'); sys.num_iters = 200;              end
    if ~isfield(sys,'noise'); sys.noise = estimate_noise(sys.y);    end
    if ~isfield(sys,'lambda'); sys.lambda = 0;
    disp('Warning: Lambda set to 0, solving the ML problem');       end
    sys.lambda = sys.lambda*sqrt(log(T)/T)*diag(sys.noise);
    if ~isfield(sys,'pen_norm'); sys.pen_norm = 'l2_sq';            end
    
    
    sys.baseline = estimate_base(sys.y,sys.noise);
    sys.y = sys.y - sys.baseline;
    sys.y(sys.y<0) = 0;
    sys.b0 = zeros(size(sys.y));
%     
%     sys.y = sys.y-repmat(min(sys.y,[],2),1,size(sys.y,2));
%     sys.b0 = 0.1*ones(size(sys.y));
    
%     sys.b0 = zeros(size(sys.y));
%     sys.b0 = repmat(min(sys.y,[],2),1,size(sys.y,2));

    if ~isfield(sys,'theta'); yt = sys.y'; 
        sys.theta = mean(aryule(yt(:,1:min(100,end)),sys.Order));   end
        % ytvec = yt(:);
        % sys.theta = aryule(ytvec(1:min(10000,end)),sys.Order);    end      
    
        sys.s0 = filter(sys.theta,1,sys.y,[],2); sys.s0(sys.s0<=0) = 1;
%     sys.s0 = ones(size(sys.y));
    
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