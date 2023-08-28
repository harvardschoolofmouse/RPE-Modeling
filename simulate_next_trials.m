% simulate_next_trials.m


if numel(ix)>1
    fast_slowMode = true;
else
    fast_slowMode = false;
end

ntrials = obj.p.weber.ntrials;
w = obj.p.weber.w;
Vh = obj.p.weber.Vh;
n = numel(Vh);
delta = obj.p.weber.delta;
n_complete = obj.p.weber.n_complete;
n_rewards = obj.p.weber.n_rewards;
t_subjective = obj.p.weber.t_subjective;
alpha = obj.params.alpha;
gamma = obj.params.gamma;
Ts = obj.p.weber.Ts;

trial = ntrials+1;
pbail = 0;



[f1,ax1] = makeStandardFigure(6,[3,2]);
C = linspecer(10); 
set(ax1(1), 'ColorOrder',C);%[blacks;C]);
set(ax1(2), 'ColorOrder',C);
set(ax1(3), 'ColorOrder',C);
title(ax1(1),'Learned')
ylabel(ax1(1),'Vh')
ylabel(ax1(3),'RPE')
ylabel(ax1(5),'Vtau (w)')
xlabel(ax1(2),'Veridical Time')
legend(ax1(1),'show','Location','Northwest')
set(f1, 'userdata', ['obj.test_RPE() CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
obj.plotLearningPerception(ax1([1,3,5]), t_subjective, Vh, delta, w, trial, [nan,nan]);


% make fast and slow clocks
% MID we will interpolate the value fxn and then trim
% we will assume Vh learned is the fast clock so that we aren't
% making data up
CS = find(t_subjective>=0,1,'first');
Vh_slow = [Vh(CS:end)';Vh(CS:end)'];
Vh_slow = [zeros(1,CS-1),reshape(Vh_slow,1,numel(Vh_slow))];
Vh_slow = Vh_slow(1:n)';
w_slow = [w(CS:end)';w(CS:end)'];
w_slow = [zeros(1,CS-1),reshape(w_slow,1,numel(w_slow))];
w_slow = w_slow(1:n)';


% make fast and slow clocks
% SLOW we will interpolate the MID value fxn and then trim
CS = find(t_subjective>=0,1,'first');
Vh_2xslow = [Vh_slow(CS:end)';Vh_slow(CS:end)'];
Vh_2xslow = [zeros(1,CS-1),reshape(Vh_2xslow,1,numel(Vh_2xslow))];
Vh_2xslow = Vh_2xslow(1:n)';
w_2xslow = [w_slow(CS:end)';w_slow(CS:end)'];
w_2xslow = [zeros(1,CS-1),reshape(w_2xslow,1,numel(w_2xslow))];
w_2xslow = w_2xslow(1:n)';

%             % FAST we will skip states and then extrapolate...
%             lasthigh = find(Vh>0,1,'last');
%             deltaVh = Vh(2:lasthigh) - Vh(1:lasthigh-1);
%             meandeltaVh = mean(deltaVh(end-5:end));
%             Vh_fast = [zeros(1,CS-1),Vh(CS:2:lasthigh), 1:(n-lasthigh).*meandeltaVh + Vh(lasthigh)];
% 
%             w_slow = [zeros(1,CS-1),reshape(w_slow,1,2*numel(w_slow))];


ixs =ix;
for ix=ixs
    T = Ts(ix);
    r = zeros(n,1); 
    if nextisrewarded
        r(T) = 1;
    end 
    disp(['rewarded? =', num2str(r(Ts(ix)))])

    if ~fast_slowMode
        [Vh,delta,w,~,~,rewarded] = tryWeberUpdate(obj,T, Vh, delta, w, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,false);
        obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh, delta, w, trial,rewarded);
    elseif fast_slowMode
        % I think we can model fast and slow clock by stretching Vh
        % and V in the weber update
        % I think this makes sense because the Vh represents the
        % subjective value estimate?
        %
        if ix <= 3 % we need relatively slow clock to pick short time
            [Vh_slow,delta,w_slow,~,~,rewarded] = tryWeberUpdate(obj,T, Vh_slow, delta, w_slow, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,false);
            obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh_slow, delta, w_slow, trial,rewarded);

        else % we need relatively fast clock to pick long time (assume the proper fxn so not making data up)
            [Vh,delta,w,~,~,rewarded] = tryWeberUpdate(obj,T, Vh, delta, w, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,false);
            obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh, delta, w, trial,rewarded);
        end
    end
end
if fast_slowMode
    title(ax1(2),'short assumed to be slow clock')
else
    title(ax1(2),['Ts=',num2str(ix)])
end