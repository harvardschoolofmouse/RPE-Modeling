% Figs 1 & 2 in Mikhael, Kim, Uchida, & Gershman.
% Written 30Sept19 by JGM.

clear; close all; clc

%% --------------------------------------------------------------------- %%

% let me start by making these params look more like my data...
T_time = 3.3;

LOI_time_ms = 400;
% cue_time_ms = 1;
% lick_time_ms = LOI_time_ms+T_time*1000;
timing_interval_ms = T_time*1000;
post_lick_interval_ms = 300;
total_post_cue_time = timing_interval_ms+post_lick_interval_ms;
total_time = LOI_time_ms+timing_interval_ms+post_lick_interval_ms;
% total_time = LOI_time_ms+cue_time_ms+timing_interval_ms+post_lick_interval_ms;

% modeled_time_ms = LOI_time_ms + lick_time_ms + 100;  % I'll try to model a lick at 3.333s, the exact reward bound

% total_time_ms = LOI_time_ms + lick_time_ms;%modeled_time_ms;  % This will give us a 400 ms pre-cue during LOI
state_width_ms = 100;                    % each veridical timing state assumed to be 100ms

% everything will be referenced to this
t_subjective = linspace(-(LOI_time_ms)/state_width_ms,...
                        total_post_cue_time/state_width_ms,...
                        total_time/state_width_ms+1)...
                        ./1000*state_width_ms;

% parameters
gamma = .9;                             % discount factor
alpha = .1;                             % learning rate
n = numel(t_subjective);    %total_time_ms/state_width_ms-1;% + LOI_time_ms/state_width_ms;     % number of states
% n_tau = (total_time_ms-LOI_time_ms)/state_width_ms;%(lick_time_ms-LOI_time_ms)/state_width_ms;% number of subjective timing states
numIter = 1000;                         % number of iterations for learning
CS = find(t_subjective==0);%LOI_time_ms/state_width_ms;        % trial onset
T_time = 3.3;%lick_time_ms/state_width_ms;%n-2;   % reward time --  this doesn't seem to make sense
weber = .15;                            % Weber fraction
T = find(t_subjective>=T_time,1,'first');
% true value
% t = 1:n_tau;%n;                                % time (ie subjective timing state), originally t by JM
t = 1:n;
r = zeros(n,1); r(T) = 1;               % reward schedule
oT = [1:CS, T+1:n];                    % times outside the trial -- also don't think we need this
V = gamma.^(T-t)'; 
% V = [V;zeros(numel(t_subjective)-numel(V),1)];
% V(1:numel(find(t_subjective<0))) = 0;
V(oT) = 0;           % true value

% [f1,ax1] = makeStandardFigure(1,[1,1]);
% plot(ax1(1),t_subjective,V)%t,xl)
% title(ax1(1),'True value landscape')
% xlabel(ax1(1),'Subjective Time')

% kernel widths
% The uncertainty kernels need to be the same length as the subjective time
% vector
% S = .1+zeros(1,n_tau);%n);                      % SD of small uncertainty kernel -- this is in PRESENCE of feedback
% L = 3+zeros(1,n_tau);%n); 
% L(lick_time_ms/state_width_ms)=0.1;% L(n-4:end) = .1;      % SD of large uncertainty kernel -- I think he sets to .1 for the end because there is feedback when reward occurs 
% web = weber*t;
S = .1+zeros(1,n);%n);                      % SD of small uncertainty kernel -- this is in PRESENCE of feedback
L = 3+zeros(1,n); L(T+1:end)=0.1;          % SD of large uncertainty kernel -- I think he sets to .1 for the end because there is feedback when reward occurs 
web = weber*(t-CS); web(1:CS)=0;

[f1,ax1] = makeStandardFigure(3,[3,1]);
plot(ax1(1),t_subjective,S)%t,xl)
plot(ax1(2),t_subjective,L)%t,xl)
plot(ax1(3),t_subjective,web)%t,xl)
title(ax1(1),'S')
title(ax1(2),'L')
title(ax1(3),'web')
xlabel(ax1(3),'Subjective Time')


% kernels = (time, kernel mean)
[xs, xl, xw] = deal(zeros(n,n));
% [xs, xl, xw] = deal(zeros(n,n-1));
% [xs, xl, xw] = deal(zeros(n_tau,n_tau-1));
% [xs, xl, xw] = deal(zeros(n,n_tau));
for y = 1:n %n_tau%n
    xs(:,y) = normpdf(t,y,S(y))';       % small kernel normpdf(x, mu, sigma)
    xl(:,y) = normpdf(t,y,L(y))';       % large kernel
    xw(:,y) = normpdf(t,y,web(y))';     % Weber's law
end
% xs = [zeros(CS-1,n_tau);xs;zeros(n-T+1,n_tau)]; xl=[zeros(CS-1,n_tau);xl;zeros(n-T+1,n_tau)]; xw=[zeros(CS-1,n_tau);xw;zeros(n-T+1,n_tau)];
% we need to append the pre-cue times to the left and top
xs = [zeros(CS-1,n_tau);xs]; xl=[zeros(CS-1,n_tau);xl]; xw=[zeros(CS-1,n_tau);xw];
xs = [zeros(CS-1,n)',xs]; xl=[zeros(CS-1,n)',xl]; xw=[zeros(CS-1,n)',xw];
% xs(:,oT)=0; xl(:,oT)=0; xw(:,oT)=0;     % leave out times outside trial
xs=xs./sum(xs); xl=xl./sum(xl); xw=xw./sum(xw);    % make prob dist's
xs(isnan(xs))=0; xl(isnan(xl))=0; xw(isnan(xw))=0; % nan's to zeros

% correction term:  [without correction | with correction]
% beta = [zeros(n,1) alpha*(exp((log(gamma))^2*(L.^2-S.^2)'/2)-1)];
% must append to L and S before putting in here
beta = [zeros(n,1) alpha*(exp((log(gamma))^2*([L(1:CS-1),L].^2-[S(1:CS-1),S].^2)'/2)-1)];

% learning with feedback
[w, Vh, delta] = deal(zeros(n,3));     	% weights, estimated value, RPE
% [w, Vh, delta] = deal(zeros(n_tau,3));     	% weights, estimated value, RPE
for e = 1:2
    for iter = 1:numIter
        for y = 1:T+1
            Vh(y,e) = w(:,e)'*xs(:,y);
            Vh(y+1,e) = w(:,e)'*xl(:,y+1);
            delta(y,e) = r(y) + gamma*Vh(y+1,e) - Vh(y,e);
            w(:,e) = w(:,e) + (alpha*delta(y,e)-beta(y,e)*w(:,e)).*xs(:,y);
            w(T+1:end,e) = 1;          % value stays high until reward
        end
    end
end

% learning without feedback
for iter = 1:numIter
    for y = 1:T+1
        Vh(y,3) = w(:,3)'*xw(:,y);
        Vh(y+1,3) = w(:,3)'*xw(:,y+1);
        delta(y,3) = r(y) + gamma*Vh(y+1,3) - Vh(y,3);
        w(:,3) = w(:,3) + alpha*delta(y,3).*xw(:,y);
        w(T+1:end,3) = r(T);         	% value stays high until reward
    end
end

%% --------------------------------------------------------------------- %%

%%
% plot uncertainty kernels
[f1,ax1] = makeStandardFigure(3,[3,1]);
plot(ax1(1),t_subjective,xl)%t,xl)
title(ax1(1),'Before Feedback')
% plot(t,xs)
plot(ax1(2),t_subjective,xs)
title(ax1(2),'After Feedback')
plot(ax1(3),t_subjective,xw)
title(ax1(3),'Weber Uncertainty')
xlabel(ax1(3),'Subjective Time')

%%
[f2,ax2] = makeStandardFigure(2,[1,2]);
fL = [3 1];                             % index of relevant value estimates
hL = [10 20 30];                        % initial location of red curves
y = 10;                                 % length of red curves along x-axis
cols=[205 92 92; 255 0 0; 150 0 0]/255; % red color schemes
feed = {'Without Feedback','With Feedback'};
for e = [2 1]
    f = fL(e);
    subplot(1,2,e) 
%     hold on
%     plot(t,V)
    plot(t_subjective,V)
%     plot(t,Vh(:,f),'k')
    plot(t_subjective,Vh(:,f),'k')
    for g = 1:3
%         redcurve = nan(1,y);
        h = hL(g);
%         redcurve(h+(0:y)) = Vh(h,f).*gamma.^(0:-1:-y);
        redcurve = Vh(h,f).*gamma.^(0:-1:-y);
        plot((hL(g)+(0:y)).*3.3/50,redcurve,'Color',cols(g,:))
    end
    % title(feed{e})
%     xlim([0 n])
    xlabel('Subjective Time')
    ylim([0 1])
    ylabel('Value')
    box off
end
legend('True Value', 'Estimated Value', 'Location','Northwest','box','off')

%%
[f3,ax3] = makeStandardFigure(4,[2,2]);
corr = {'Without Correction','With Correction'};
for e = 1:2
    subplot(2,2,e)
%     hold on
    t_subjective
%     plot(t,V)
    plot(t_subjective,V)
%     plot(t,Vh(:,e),'k')
    plot(t_subjective,Vh(:,e),'k')
    ylabel('Value')
    ylim([min(min([Vh V])) max(max([Vh V]))])
    title(corr{e})
    
    subplot(2,2,e+2)
%     plot(t,delta(:,e),'k')
    plot(t_subjective,delta(:,e),'k')
    ylabel('RPE')
    ylim([min(min(delta(:,1:2))) max(max(delta(:,1:2)))])
end

subplot(2,2,1)
legend('True Value','Estimated Value','Location','Northwest','box','off')

for e = 1:4
    subplot(2,2,e)
    xlabel('Time')
%     xlim([0 n])
    box off
end