function sys_fade = FADE( sys )
sys = default_opts_fade(sys);
[p,T] = size(sys.y);

b = sys.b0;
s = sys.s0;
% objective = zeros(1,sys.num_iters);
ds = zeros(1,sys.num_iters);
iter = 1;
ds(1) = 1;
    while iter < sys.num_iters && ( ds(iter)>0.005 || iter<sys.min_iters)
        [Lp,Ln] = calc_L(sys,s,b);
        Gp = calc_penalty(sys,s);
%         X = filter(1,sys.theta,s,[],2);
%         bnew = min(X,[],2);
        
%         bnew = max(sum(sys.y-X,2)/T,0);

%         bnew = sum(sys.y,2)./sum(X+b,2);
        snew = Ln./(Lp+sys.lambda*Gp).*s;


%         objective(iter) = norm(diag(1./sys.noise.^2)*(sys.y-X),'fro')^2 + sys.lambda*sum(s(:));
        ds(iter+1) = max(abs((sum(s,2)-sum(snew,2))./sum(s,2)));
        
        s = snew;
%         b = repmat(bnew,1,T);
        iter = iter + 1;
    end
X = filter(1,sys.theta,s,[],2);

sys_fade = sys;
% sys_fade.b = b;
sys_fade.spikes = s;
sys_fade.smoothed_traces = X + sys.baseline;
sys_fade.theta = sys.theta;
% sys_fade.objective = objective(1:iter-1);
sys_fade.ds = ds(1:iter-1);

end

function [Lp,Ln] = calc_L(sys,s,b)
Ln = 2*fliplr(filter(1,sys.theta,fliplr(sys.y),[],2));
X = filter(1,sys.theta,s,[],2);
Lp = 2*fliplr(filter(1,sys.theta,fliplr(X+b),[],2));
end

function Gp = calc_penalty(sys,s)
    switch sys.pen_norm
        case 'l1'
            Gp = ones(size(s));
        case 'l2_sq'
            Gp = 2*s;
        case 'l2'
            T = size(s,2);
            ss = repmat(sqrt(sum(s.^2,2)),1,T);
            Gp = 2*s./ss;
    end

end