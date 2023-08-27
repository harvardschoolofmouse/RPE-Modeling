% behavioral_timing_value_landscape.m
%   Created 12/24/21 - Allison Hamilos (ahamilos [at] mit [dot] edu)
%   Based on original code from John Mikhael (Mikhael, Kim, Uchida, &
%   Gershman.)

clear; close all; clc

%% --------------------------------------------------------------------- %%
debug=false;

T_time = 3.3;           % target time for optimal timer, this is when lick occurs. Reward time in Mikhael model

% intervals -- use to construct the subjective timespace (statespace)
LOI_time_ms = 400;                      % lamp-off interval
timing_interval_ms = T_time*1000;       % timing interval from cue to lick (not inclusive of 0)
post_lick_interval_ms = 300;            % time to look post-lick
total_post_cue_time = timing_interval_ms+post_lick_interval_ms;     % total states past the cue (non inclusive of 0)
total_time = LOI_time_ms+timing_interval_ms+post_lick_interval_ms;  % total time (not inclusive of 0)
state_width_ms = 80;                    % each subjective timing state assumed to include this many ms of veridical time

% everything will be referenced to the subjective state space, e.g.,
%   [-0.4s:0.1s:3.3s+post-lick time]
t_subjective = linspace(-(LOI_time_ms)/state_width_ms,...
                        total_post_cue_time/state_width_ms,...
                        total_time/state_width_ms+1)...
                        ./1000*state_width_ms;

% TD learning parameters
gamma = .9;                             % discount factor
alpha = .1;                             % learning rate
n = numel(t_subjective);                % number of subjective timing states
numIter = 1000;                         % number of iterations for learning
CS = find(t_subjective==0);             % cue position (by def'n, 0 s)
T = find(t_subjective>=T_time,1,'first'); % lick state position
weber = .15;                            % Weber fraction
% t_indexed = [zeros(1,CS),1:sum(t_subjective>0)];    % we'll use this to deal with value learning of
                                        % the timing states. There's no
                                        % info about the state until the
                                        % cue comes along

% true value landscape         
t = 1:n;                                % time (ie subjective timing state), originally t by JM
r = zeros(n,1); r(T) = 1;               % reward schedule (i.e., position of lick)
oT = [1:CS, T+1:n];                     % state positions outside the trial (cue and before, post-lick)
V = gamma.^(T-t)'; V(oT) = 0;           % assetion: true value landscape via exponential TD model
if debug
    [f1,ax1] = makeStandardFigure(1,[1,1]);
    plot(ax1(1),t_subjective,V)%t,xl)
    title(ax1(1),'True value landscape')
    xlabel(ax1(1),'Subjective Time')
end

% state space uncertainty kernel widths
%
%   under Mikhael model, ramping will result from resolved state
%   uncertainty. the animal will impose its own resolution because
%   it must make a decision
%
S = .1+zeros(1,n);                      % SD of small uncertainty kernel -- this is in PRESENCE of feedback
L = 3+zeros(1,n); L(T-3:end)=.1;%.1;       % SD of large uncertainty kernel -- drops to 0.1 by Mikhael convention after reward feedback 
web = weber*(t-CS); web(1:CS)=0;        % we assert there is no weber uncertainty before the cue, or rather uncertainty is max here

if debug
    [f1,ax1] = makeStandardFigure(3,[3,1]);
    plot(ax1(1),t_subjective,S)%t,xl)
    plot(ax1(2),t_subjective,L)%t,xl)
    plot(ax1(3),t_subjective,web)%t,xl)
    title(ax1(1),'S')
    title(ax1(2),'L')
    title(ax1(3),'web')
    xlabel(ax1(3),'Subjective Time')
end

% construction of state space uncertainty kernels
%       kernels = (probability of being being in veridical t, subjective state)
%       normpdf(veridical time t, average subjective time tau when in veridical state t, kernel width at veridical time t)
%  These are the probability of being in state t|subjective state tau. 
%       [tau x t]
%   xs = p(t|tau, sigma=small=0.1)
%   xl = p(t|tau, sigma=small=0.1)
%   xs = p(t|tau, sigma=small=0.1)
%
[xs, xl, xw] = deal(zeros(n,n));
if debug
    [f1,ax1] = makeStandardFigure(3,[3,1]);
    title(ax1(1),'p(t|tau, \sigma_s')
    title(ax1(2),'p(t|tau, \sigma_l * \sigma_{weber})')
    title(ax1(3),'p(t|tau, \sigma_{weber})')
    xlabel(ax1(3),'Veridical Time (s)')
end
for y = 1:n 
    xs(:,y) = normpdf(t,y,S(y))';           % small kernel, 
    xl(:,y) = normpdf(t,y,L(y))';           % large kernel
    xw(:,y) = normpdf(t,y,web(y))';     % Weber's law

end
% xl = xl.*xw;
% xs = xs;
xs(:,oT)=0; xl(:,oT)=0; xw(:,oT)=0;     % leave out times outside trial
xs=xs./sum(xs); xl=xl./sum(xl); xw=xw./sum(xw);    % make prob dist's
xs(isnan(xs))=0; xl(isnan(xl))=0; xw(isnan(xw))=0; % nan's to zeros

if debug
    plot(ax1(1),t_subjective,xs)
    plot(ax1(2),t_subjective,xl)
    plot(ax1(3),t_subjective,xw)
end 

% correction term:  [without correction | with correction]
% this term updates the estimate of the value of veridical state t (Vh_t)
% based on the current estimate of Vh_t + learning rate applied to the RPE
% at subjective state tau (times probability of being in veridical t given
% tau) - corrected (beta) estimate of value at tau times prob of being in t
% given tau
beta = [zeros(n,1), alpha*(exp((log(gamma))^2*(L.^2-S.^2)'/2)-1)];
% the correction will be applied anytime sensory feedback is provided

% we assume that the animal is constantly resolving its own state
% uncertainty by virtue of the fact it must act upon its estimate of time

% learning with feedback
[w,... the learned relevance of each state (the weights of each subjective timing state)
    Vh,... the estimate of the value function at subjective timing state tau
    delta,... the rpe at subjective timing state tau
    ] = deal(zeros(n,3));              % weights, estimated value, RPE
[f1,ax1] = makeStandardFigure(2,[2,1]);
C = linspecer(10); 
set(ax1(1), 'ColorOrder',C);%[blacks;C]);
set(ax1(2), 'ColorOrder',C);
title(ax1(1),['V_h, beta=', num2str(beta(1,1))])
title(ax1(2),'RPE')
xlabel(ax1(2),'Subjective Time')
%
% NB! w = Vh_tau, the subjective estimate of value. it's NOT the weight of
% the state
%

for e = 1:1 % beta=0, no feedback
    for iter = 1:numIter % for each iteration of learning...
        for y = 1:T+1
            Vh(y,e) = w(:,e)'*xs(:,y); % Vh=Vt, w=Vtau, xs=p(t|tau, sigma=s) -- assuming we got feedback
            Vh(y+1,e) = w(:,e)'*xl(:,y+1); % Vh=Vt+1, w=Vtau+1, xs=p(t+1|tau+1, sigma=l)
            delta(y,e) = r(y) + gamma*Vh(y+1,e) - Vh(y,e);
            w(:,e) = w(:,e) + (alpha*delta(y,e)-beta(y,e)*w(:,e)).*xs(:,y); 
            w(T+1:end,e) = 1;          % value stays high until reward
        end
%         if debug
        if ismember(iter,(numIter/10).*[1:10])
            plot(ax1(1),t_subjective,Vh(:,e), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
            plot(ax1(2),t_subjective,delta(:,e), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
%         end
        end
    end
end
legend(ax1(1),'show')
[f1,ax1] = makeStandardFigure(2,[2,1]);
title(ax1(1),['V_h, beta=', num2str(beta(1,2))])
title(ax1(2),'RPE')
C = linspecer(10); 
set(ax1(1), 'ColorOrder',C);%[blacks;C]);
set(ax1(2), 'ColorOrder',C);
xlabel(ax1(2),'Subjective Time')
for e = 2:2 % beta=beta, with feedback
    for iter = 1:numIter % for each iteration of learning...
        for y = 1:T+1
            Vh(y,e) = w(:,e)'*xs(:,y);
            Vh(y+1,e) = w(:,e)'*xl(:,y+1);
            delta(y,e) = r(y) + gamma*Vh(y+1,e) - Vh(y,e);
            w(:,e) = w(:,e) + (alpha*delta(y,e)-beta(y,e)*w(:,e)).*xs(:,y);
            w(T+1:end,e) = 1;          % value stays high until reward
        end
        if ismember(iter,(numIter/10).*[1:10])
            plot(ax1(1),t_subjective,Vh(:,e), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
            plot(ax1(2),t_subjective,delta(:,e), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
        end
    end
end
legend(ax1(1),'show')

% learning without feedback
[f1,ax1] = makeStandardFigure(2,[2,1]);
title(ax1(1),['V_h, no feedback (weber uncertainty)'])
title(ax1(2),'RPE')
xlabel(ax1(2),'Subjective Time')
C = linspecer(10); 
set(ax1(1), 'ColorOrder',C);%[blacks;C]);
set(ax1(2), 'ColorOrder',C);
for iter = 1:numIter
    for y = 1:T+1
        Vh(y,3) = w(:,3)'*xw(:,y);
        Vh(y+1,3) = w(:,3)'*xw(:,y+1);
        delta(y,3) = r(y) + gamma*Vh(y+1,3) - Vh(y,3);
        w(:,3) = w(:,3) + alpha*delta(y,3).*xw(:,y);
        w(T+1:end,3) = r(T);         	% value stays high until reward
    end
    if ismember(iter,(numIter/10).*[1:10])
        plot(ax1(1),t_subjective,Vh(:,3), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
        plot(ax1(2),t_subjective,delta(:,3), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
    end
end
legend(ax1(1),'show')

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
hL = round([0.25*3300/state_width_ms 0.5*3300/state_width_ms 0.75*3300/state_width_ms]);                        % initial location of red curves
y = round(0.25*3300/state_width_ms);                                 % length of red curves along x-axis
cols=[205 92 92; 255 0 0; 150 0 0]/255; % red color schemes
feed = {['State uncertainty but no feedback (beta=', num2str(beta(1,1)),')'],['State uncertainty with continual feedback beta=(', num2str(beta(1,2)) ')']};
for e = [2 1]
    f = fL(e);
    subplot(1,2,e) 
    plot(t_subjective,V, 'linewidth', 3)
    plot(t_subjective,Vh(:,f),'k', 'linewidth', 3)
    for g = 1:3
        redcurve = nan(1,n);
%         redcurve = nan(1,y);
        h = hL(g);
        redcurve(h+(0:y)) = Vh(h,f).*gamma.^(0:-1:-y);
%         redcurve = Vh(h,f).*gamma.^(0:-1:-y);
%         plot((hL(g)+(0:y)).*3.3/50,redcurve,'Color',cols(g,:))
        plot(t_subjective,redcurve,'Color',cols(g,:), 'linewidth', 3)
    end
    title(feed{e})
%     xlim([0 n])
    xlabel('Subjective Time')
    ylim([0 1])
    ylabel('Value')
    box off
end
legend('True Value', 'Estimated Value', 'Location','Northwest','box','off')

%%
[f3,ax3] = makeStandardFigure(6,[2,3]);
corr = {'Without Correction','With Correction','Weber'};
for e = 1:2
    plot(ax3(e),t_subjective,V, 'linewidth', 3)
    plot(ax3(e),t_subjective,Vh(:,e),'k', 'linewidth', 3)
    ylabel(ax3(e),'Value')
    ylim(ax3(e),[min(min([Vh V])) max(max([Vh V]))])
    title(ax3(e),corr{e})
    
    
    plot(ax3(e+3),t_subjective,delta(:,e),'k', 'linewidth', 3)
    ylabel(ax3(e+3),'RPE')
    ylim(ax3(e+3),[min(min(delta(:,1:2))) max(max(delta(:,1:2)))])
end

plot(ax3(3),t_subjective,V, 'linewidth', 3)
plot(ax3(3),t_subjective,Vh(:,3),'k', 'linewidth', 3)
ylabel(ax3(3),'Value')
ylim(ax3(3),[min(min([Vh V])) max(max([Vh V]))])
title(ax3(3),corr{3})

plot(ax3(6),t_subjective,delta(:,3),'k', 'linewidth', 3)
ylabel(ax3(6),'RPE')
ylim(ax3(6),[min(min(delta(:,3))) max(max(delta(:,3)))])

legend(ax3(1),'True Value','Estimated Value','Location','Northwest','box','off')

for e = 1:6
    xlabel(ax3(e),'Time')
    box off
end