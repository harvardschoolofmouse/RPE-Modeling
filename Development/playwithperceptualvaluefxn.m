%
% initialize standard params
%   a major difference here is that we have feedback coming at different
%   times and reward comes after that...
%   with feedback, we must exit the task...
%
T_time = [0.6, 1.05, 1.26, 1.74, 1.95, 2.4];           % target time for optimal timer, this is when lick occurs. Reward time in Mikhael model
p_correct = [0.98, 0.94, 0.82, 0.71, 0.85, 0.94];

% intervals -- use to construct the subjective timespace (statespace)
precue_time_ms = 400;                   % time before cue
timing_interval_ms = T_time.*1000;      % timing intervals from cue to stop
post_stop_interval_ms = 600;            % time to look after stop cue
total_post_cue_time = timing_interval_ms+post_stop_interval_ms;     % total states past the cue (non inclusive of 0)
total_time = precue_time_ms+timing_interval_ms+post_stop_interval_ms;  % total time for each interval (not inclusive of 0)
state_width_ms = 80;                    % each subjective timing state assumed to include this many ms of veridical time

% everything will be referenced to the subjective state space, e.g.,
%   we will have a separate timeline for each case
t_subjective = nan(6, round(total_time(end)/state_width_ms+1));
for ii = 1:numel(T_time)
    t_subjective(ii,:) = linspace(-precue_time_ms/state_width_ms,...
                        total_post_cue_time(end)/state_width_ms,...
                        round(total_time(end)/state_width_ms+1))...
                        ./1000*state_width_ms;   
    t_subjective(ii,t_subjective(ii,:)>total_post_cue_time(ii)/1000) = nan;
end

% TD learning parameters
gamma = .9;                             % discount factor
alpha = .1;                             % learning rate
n = numel(t_subjective(1,:));                % number of subjective timing states
numIter = 1000;                         % number of iterations for learning
CS = find(t_subjective==0);             % cue position (by def'n, 0 s)
T = find(t_subjective>=T_time,1,'first'); % lick state position
weber = .15;                            % Weber fraction


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