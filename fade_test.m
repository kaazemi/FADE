clear;clc; close all;
p = 100; T = 10000; s = 100;
w = zeros(p,T);
for i=1:p
w(i,randsample(T,s)) = 1;
end
theta = [1 -.95];
b = 2;
X = filter(1,theta,w,[],2);
sys.y = X+0.1*randn(size(w))+b;
%%
% sys.Order = 2;
sys.lambda = 1000;
% sys.pen_norm = 'l1';
sys.p_norm = 1;
sys.q_norm = 1;
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
f3 = subplot(3,1,3); stem(sys_fade.spikes(index,:),'marker','none'); title('Recovered')
linkaxes([f1 f2 f3],'x')
%%%%%%%
% figure(2);
% plot(sys_fcss.objective(2:end));
figure;
plot(sys_fade.ds(2:end));
% sys_fcss.b
%%%%%%