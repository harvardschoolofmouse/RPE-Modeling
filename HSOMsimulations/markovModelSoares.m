trials = 1000;

ms_per_state = 100;
fb_idx_wrt_timing = [find(1:(2400/ms_per_state) >= 600/ms_per_state,1,'first'),...
                    find(1:(2400/ms_per_state) >= 1050/ms_per_state,1,'first'),...
                    find(1:(2400/ms_per_state) >= 1260/ms_per_state,1,'first'),...
                    find(1:(2400/ms_per_state) >= 1740/ms_per_state,1,'first'),...
                    find(1:(2400/ms_per_state) >= 1950/ms_per_state,1,'first'),...
                    find(1:(2400/ms_per_state) >= 2400/ms_per_state,1,'first')];
p_correct = [0.98, 0.94, 0.82, 0.71, 0.85, 0.94];

s = 0.1;
l = 3;
weber = 0.15;

% start by defining states, assume a perfect timer
%pre-cue states
states_precue = 5;
p_txn_precue = ones(1,5);
v_precue = zeros(1,5);
delta_precue = zeros(1,5);
% cue state
states_cue = nan;
p_txn_cue = 1;
v_cue = 0;
delta_cue = 0;
% timing states
states_timing = 2400/ms_per_state;
p_bail = ones(1,states_timing).*(1/(2*states_timing));
p_txn_timing = ones(1, 2400/ms_per_state); p_txn_timing(fb_idx_wrt_timing) = 1/6; p_txn_timing = p_txn_timing-p_bail;
v_timing = zeros(1,states_timing);
delta_cue = zeros(1,states_timing);
% feedback states
states_fb = 300/ms_per_state;
p_txn_fb = ones(1,states_fb);
% reward/fail states
states_reward = 1;
p_txn_reward = 1;
states_fail = 1;
p_txn_fail = 1;


% all states and their uncertainty kernels
allstates_precue_idx = 1:5;
allstates_cue_idx = 6;
allstates_timing_idx = 7:7+states_timing;


allstates = {
    'precue 1';...
    'precue 2';...
    'precue 3';...
    'precue 4';...
    'precue 5';...
    'cue'};
for ii = 1:states_timing
    allstates{end+1} = ['timing ',num2str(ii)];
end

for fb = 1:6
    allstates_fb_idx(fb, 1:numel(numel(allstates)+1:numel(allstates)+states_fb)) = numel(allstates)+1:numel(allstates)+states_fb;
    for ii = 1:states_fb
        allstates{end+1} = ['fb',num2str(fb), ' ', num2str(ii)];
    end
end

allstates{end+1} = 'reward';
allstates_reward_idx = numel(allstates);
allstates{end+1} = 'fail';
allstates_fail_idx = numel(allstates);


visits = zeros(numel(allstates),1);
% run simulation
state = 'precue';
n_in_state = 1;
t=0;
ntrials = 1000;
nattempts = zeros(6,1);
ncompletions = zeros(6,1);
nrewards = zeros(6,1);
nbails = 0;
[f,ax] = makeStandardFigure();
while t < ntrials+1
    switch state
        case 'precue'
            if n_in_state==1
                if ismember(t, 1000)
                    cla(ax)
                    bar(ax,nattempts,'k', 'displayname','attempts'), bar(ax,ncompletions,'displayname','completed'), bar(ax,nrewards,'g','displayname','rewarded')
                    legend('show')
                    xlabel('stop cue time #')
                    title(['trial#' num2str(t)])
                    ylabel('# trials')
                end
                t = t+1;
                trial_type = randi(6);
                nattempts(trial_type) = nattempts(trial_type) + 1;
            end
            if n_in_state < states_precue
                if rand<p_txn_precue(n_in_state)
                    visits(allstates_precue_idx(n_in_state)) = visits(allstates_precue_idx(n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'precue';
                end
            elseif rand<p_txn_precue(n_in_state)
                visits(allstates_precue_idx(n_in_state)) = visits(allstates_precue_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'cue';
            else
                error('');
            end
        case 'cue'
            if n_in_state < states_cue
                if rand<p_txn_cue(n_in_state)
                    visits(allstates_cue_idx(n_in_state)) = visits(allstates_cue_idx(n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'cue';
                end
            elseif rand<p_txn_cue(n_in_state)
                visits(allstates_cue_idx(n_in_state)) = visits(allstates_cue_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'timing';
            end
        case 'timing'
            if n_in_state==fb_idx_wrt_timing(1) && trial_type == 1
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb1';
                continue
            elseif n_in_state==fb_idx_wrt_timing(2) && trial_type == 2
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb2';
                continue
            elseif n_in_state==fb_idx_wrt_timing(3) && trial_type == 3
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb3';
                continue
            elseif n_in_state==fb_idx_wrt_timing(4) && trial_type == 4
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb4';
                continue
            elseif n_in_state==fb_idx_wrt_timing(5) && trial_type == 5
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb5';
                continue
            elseif n_in_state==fb_idx_wrt_timing(6) && trial_type == 6
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                n_in_state = 1;
                state = 'fb6';
                continue
            elseif n_in_state < states_timing
                if rand<p_bail(n_in_state)
                    visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                    n_in_state = 1;
                    nbails = nbails+1;
                    state = 'precue';
                elseif rand<p_txn_timing(n_in_state)
                    visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'timing';
                end
            elseif rand<p_txn_timing(n_in_state)
                visits(allstates_timing_idx(n_in_state)) = visits(allstates_timing_idx(n_in_state))+1;
                warning('got to cat 6 without a flip...')
                n_in_state = 1;
                state = 'fb6';
            end
        case 'fb1'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(1,n_in_state)) = visits(allstates_fb_idx(1,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb1';
                end
            elseif rand<p_correct(1)
                visits(allstates_fb_idx(1,n_in_state)) = visits(allstates_fb_idx(1,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(1,n_in_state)) = visits(allstates_fb_idx(1,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'fb2'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(2,n_in_state)) = visits(allstates_fb_idx(2,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb2';
                end
            elseif rand<p_correct(2)
                visits(allstates_fb_idx(2,n_in_state)) = visits(allstates_fb_idx(2,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(2,n_in_state)) = visits(allstates_fb_idx(2,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'fb3'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(3,n_in_state)) = visits(allstates_fb_idx(3,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb3';
                end
            elseif rand<p_correct(3)
                visits(allstates_fb_idx(3,n_in_state)) = visits(allstates_fb_idx(3,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(3,n_in_state)) = visits(allstates_fb_idx(3,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'fb4'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(4,n_in_state)) = visits(allstates_fb_idx(4,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb4';
                end
            elseif rand<p_correct(4)
                visits(allstates_fb_idx(4,n_in_state)) = visits(allstates_fb_idx(4,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(4,n_in_state)) = visits(allstates_fb_idx(4,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'fb5'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(5,n_in_state)) = visits(allstates_fb_idx(5,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb5';
                end
            elseif rand<p_correct(5)
                visits(allstates_fb_idx(5,n_in_state)) = visits(allstates_fb_idx(5,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(5,n_in_state)) = visits(allstates_fb_idx(5,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'fb6'
            if n_in_state < states_fb
                if rand<p_txn_fb(n_in_state)
                    visits(allstates_fb_idx(6,n_in_state)) = visits(allstates_fb_idx(6,n_in_state))+1;
                    n_in_state = n_in_state + 1;
                    state = 'fb6';
                end
            elseif rand<p_correct(6)
                visits(allstates_fb_idx(6,n_in_state)) = visits(allstates_fb_idx(6,n_in_state))+1;
                n_in_state = 1;
                state = 'reward';
            else
                visits(allstates_fb_idx(6,n_in_state)) = visits(allstates_fb_idx(6,n_in_state))+1;
                n_in_state = 1;
                state = 'fail';
            end
        case 'reward'
            if n_in_state < states_reward
                if rand<p_txn_reward(n_in_state)
                    visits(allstates_reward_idx-1+n_in_state) = visits(allstates_reward_idx-1+n_in_state)+1;
                    n_in_state = n_in_state + 1;
                    state = 'reward';
                end
            else
                visits(allstates_reward_idx-1+n_in_state) = visits(allstates_reward_idx-1+n_in_state)+1;
                n_in_state = 1;
                ncompletions(trial_type) = ncompletions(trial_type)+1;
                nrewards(trial_type) = nrewards(trial_type)+1;
                state = 'precue';
            end
        case 'fail'
            if n_in_state < states_fail
                if rand<p_txn_fail(n_in_state)
                    visits(allstates_fail_idx-1+n_in_state) = visits(allstates_fail_idx-1+n_in_state)+1;
                    n_in_state = n_in_state + 1;
                    state = 'fail';
                end
            else
                visits(allstates_fail_idx-1+n_in_state) = visits(allstates_fail_idx-1+n_in_state)+1;
                n_in_state = 1;
                ncompletions(trial_type) = ncompletions(trial_type)+1;
                state = 'precue';
            end
        otherwise
            error('')
    end
end


% allstates_fb1_idx = 7+allstates_timing_idx(fb_idx_wrt_timing(1):7+
% obj.b.params.xs(:,y) = normpdf(obj.b.params.t,y,obj.b.params.S(y))'; 





    