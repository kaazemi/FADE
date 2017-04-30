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
        Pp = calc_penalty(sys,s);
%         X = filter(1,sys.theta,s,[],2);
%         bnew = min(X,[],2);
        
%         bnew = max(sum(sys.y-X,2)/T,0);

%         bnew = sum(sys.y,2)./sum(X+b,2);
        snew = Ln./(Lp+sys.lambda*Pp).*s;


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

function Pp = calc_penalty(sys,s)
    A1 = sum(sum(s.^sys.p_norm,2).^(sys.q_norm/sys.p_norm))^(1/sys.q_norm-1);
    A2 = repmat(sum(s.^sys.p_norm,2).^(sys.q_norm/sys.p_norm-1),1,sys.T);
    A3 = s.^(sys.p_norm-1);
    Pp = A1*A2.*A3;
%     switch sys.pen_norm
%         case 'l1'
%             Pp = ones(size(s));
%         case 'l2_sq'
%             Pp = 2*s;
%         case 'l2'
%             T = size(s,2);
%             ss = repmat(sqrt(sum(s.^2,2)),1,T);
%             Pp = 2*s./ss;
%     end

end










