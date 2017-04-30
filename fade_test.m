clear;clc; close all;
%% Generate data
p = 100;   % 100 ROI's
T = 10000; % number of frames
s = 100;   % sparsity 
w = zeros(p,T); %spikes
for i=1:p
w(i,randsample(T,s)) = 1;
end
theta = [1 -.95]; % AR parameters
b = 2; % baseline
X = filter(1,theta,w,[],2); %clean calcium traces
sys.y = X+0.1*randn(size(w))+b; %Noisy calcium traces
%%
sys.lambda = 1000; %regularization parameter
sys.p_norm = 1;    % p norm in ell_pq penalty on spikes
sys.q_norm = 1;    % q norm in ell_pq penalty on spikes
message = sprintf('Penalizing ell_%d,%d norm of spikes',sys.p_norm,sys.q_norm);
disp(message)

tic
sys_fade = FADE( sys );
toc
%% 
index = 1;
figure(1);
f1 = subplot(3,1,1); plot(sys.y(index,:)); title('Observed')
f2 = subplot(3,1,2); plot(sys_fade.smoothed_traces(index,:)); title('Ground-Truth')
hold on; stem(w(index,:),'marker','none');
f3 = subplot(3,1,3); stem(sys_fade.spikes(index,:),'marker','none'); title('Recovered Spikes: No post processing done')
linkaxes([f1 f2 f3],'x')
%%%%%%%
% figure(2);
% plot(sys_fcss.objective(2:end));
figure;
plot(100*sys_fade.ds(2:end)); title('Relative Changes in estimated spikes')
ylabel('Percent'), xlabel('Iteration')
% sys_fcss.b
%%%%%%