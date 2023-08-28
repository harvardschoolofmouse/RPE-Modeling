classdef CLASS_value_landscaping_obj < handle
	% 
	% 	Created 	12/25/21	ahamilos
	% 	Modified 	08/27/23	ahamilos
	% 
	% 	Based on simulations originally built and reported by John Mikhael et al, 2019
	% 
	properties
		params
        b
        p
    end

	
	methods
		%-------------------------------------------------------
		%		Methods: Initialization
		%-------------------------------------------------------
		function obj = CLASS_value_landscaping_obj()
            obj.params.runID = randi(10000);
            %
            % initialize standard params
            %
			T_time = 3.3;           % target time for optimal timer, this is when lick occurs. Reward time in Mikhael model

            % intervals -- use to construct the subjective timespace (statespace)
            LOI_time_ms = 400;                      % lamp-off interval
            timing_interval_ms = T_time*1000;       % timing interval from cue to lick (not inclusive of 0)
            post_lick_interval_ms = 1100;%300;            % time to look post-lick
            warning('change here')
            total_post_cue_time = timing_interval_ms+post_lick_interval_ms;     % total states past the cue (non inclusive of 0)
            total_time = LOI_time_ms+timing_interval_ms+post_lick_interval_ms;  % total time (not inclusive of 0)
            state_width_ms = 80;                    % each subjective timing state assumed to include this many ms of veridical time

            % everything will be referenced to the subjective state space, e.g.,
            %   [-0.4s:0.1s:3.3s+post-lick time]
            t_subjective = linspace(-(LOI_time_ms)/state_width_ms,... % this is pre-cue time
                                    total_post_cue_time/state_width_ms,... % this is total time post-cue
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

            obj.params.gamma = gamma;
            obj.params.alpha = alpha;
            obj.params.numIter = numIter;
            obj.params.weber = weber;
            

            obj.b.params.n = n;
            obj.b.params.CS = CS;
            obj.b.params.T = T;
            obj.b.params.T_time = T_time;
            obj.b.params.LOI_time_ms = LOI_time_ms;
            obj.b.params.timing_interval_ms = timing_interval_ms;
            obj.b.params.post_lick_interval_ms = post_lick_interval_ms;
            obj.b.params.total_post_cue_time = total_post_cue_time;
            obj.b.params.total_time = total_time;
            obj.b.params.state_width_ms = state_width_ms;
            obj.b.params.t_subjective = t_subjective;
            
            % set params
            obj.makeBehavioralValueLandscape();
            obj.defineUncertaintyKernels();
            obj.defineBeta();
            
            % do learning with Mikhael et al model
            obj.learning(1);
            obj.learning(2);
            obj.learning(3);
            
            % simulate difference pacemaker rates and the corresponding
            % value trajectories
            obj.simulateBehavior();
            obj.simBehaviorValueTrajectory()
        end

        function makeBehavioralValueLandscape(obj,Plot)
            if nargin < 2
                Plot=false;
            end
            % true value landscape         
            obj.b.params.t = 1:obj.b.params.n;      % time (ie subjective timing state), originally t by JM
            warning('changed this. mouse doesnt know outcome till after lick')
            obj.b.params.r = zeros(obj.b.params.n,1); obj.b.params.r(obj.b.params.T+1) = 1;   % reward schedule (i.e., position of lick)
            obj.b.params.oT = [1:obj.b.params.CS, obj.b.params.T+1:obj.b.params.n];   % state positions outside the trial (cue and before, post-lick)
            obj.b.params.V = obj.params.gamma.^(obj.b.params.T-obj.b.params.t)'; 
%             obj.b.params.V(obj.b.params.oT) = 0;           % assetion: true value landscape via exponential TD model
            warning('changed this. trying to fix discontinuity')
            obj.b.params.V(1:obj.b.params.CS) = 0; % assertion -- value after the lick is not defined, so continues to increase
            if Plot
                [f1,ax1] = makeStandardFigure(1,[1,1]);
                plot(ax1(1),obj.b.params.t_subjective,obj.b.params.V)
                title(ax1(1),'True value landscape')
                xlabel(ax1(1),'Veridical Time (s)')
                ylabel(ax1(1),'V_t')
                set(f1, 'userdata', ['makeBehavioralValueLandscape(obj,Plot) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            end
        end
        function makePerceptualValueLandscape(obj,Plot)
            if nargin < 2
                Plot=false;
            end
            %
            % initialize standard params
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
            
            % true value landscape         
            obj.p.params.n = size(t_subjective,2);
            obj.p.params.t = 1:obj.b.params.n;      % time (ie subjective timing state), originally t by JM
            obj.p.params.r = zeros(obj.b.params.n,1); obj.b.params.r(obj.b.params.T) = 1;               % reward schedule (i.e., position of lick)
            obj.p.params.oT = [1:obj.b.params.CS, obj.b.params.T+1:obj.b.params.n];                     % state positions outside the trial (cue and before, post-lick)
            obj.p.params.V = obj.params.gamma.^(obj.b.params.T-obj.b.params.t)'; obj.b.params.V(obj.b.params.oT) = 0;           % assetion: true value landscape via exponential TD model
            if Plot
                [f1,ax1] = makeStandardFigure(1,[1,1]);
                plot(ax1(1),obj.b.params.t_subjective,obj.b.params.V)
                title(ax1(1),'True value landscape')
                xlabel(ax1(1),'Veridical Time (s)')
                ylabel(ax1(1),'V_t')
                set(f1, 'userdata', ['makeBehavioralValueLandscape(obj,Plot) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            end
        end
        function sts = smooth(obj, ts, OnOff, method)
			% 
			% OnOff = 0: no smoothing
			% OnOff > 0: the kernel is OnOff
			% 
			if nargin < 4
				method = 'gausssmooth';
			end
			if nargin < 3 || OnOff < 0
				OnOff = obj.Plot.smooth_kernel;
            end
            
            if isempty(ts), sts = []; return, end

			if strcmp(method, 'gausssmooth')
				if OnOff
					sts = gausssmooth(ts, round(OnOff), 'gauss');
				else
					sts = ts;
				end
			else
				if OnOff
					sts = smooth(ts, round(OnOff), 'moving');
				else
					sts = ts;
				end
			end
		end
        function defineUncertaintyKernelsWidths(obj,Plot,s,l,w)
            if nargin < 2
                Plot=true;
            end
            if nargin < 3
                s=0.1;
                l=3;
                w=obj.params.weber;
            end
            obj.params.kernels.s = s;
            obj.params.kernels.l = l;
            obj.params.kernels.w = w;
            % state space uncertainty kernel widths
            %
            %   under Mikhael model, ramping will result from resolved state
            %   uncertainty. the animal will impose its own resolution because
            %   it must make a decision
            %
            % behavioral mode
            obj.b.params.S = s+zeros(1,obj.b.params.n);    % SD of small uncertainty kernel -- this is in PRESENCE of feedback
            obj.b.params.L = l+zeros(1,obj.b.params.n); 
%             obj.b.params.L(obj.b.params.T-1:end)=.1;% we will assume animal knows it is licking just before the reward state
            obj.b.params.web = w*(obj.b.params.t-obj.b.params.CS); 
            obj.b.params.web(1:obj.b.params.CS)=0;        % we assert there is no weber uncertainty before the cue, or rather uncertainty is max here
            % perceptual mode
            %
            if Plot
                [f1,ax1] = makeStandardFigure(3,[3,1]);
                plot(ax1(1),obj.b.params.t_subjective,obj.b.params.S)
                plot(ax1(2),obj.b.params.t_subjective,obj.b.params.L)
                plot(ax1(3),obj.b.params.t_subjective,obj.b.params.web)
                title(ax1(1),'S')
                title(ax1(2),'L')
                title(ax1(3),'web')
                xlabel(ax1(3),'Objective Time')
                set(f1, 'userdata', ['defineUncertaintyKernelsWidths(obj,Plot,s,l,w) s=', num2str(s), ' l=', num2str(l), ' w=', num2str(w), ' CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            end
        end
        function defineUncertaintyKernels(obj,Plot)
            if nargin < 2
                Plot=true;
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
            % behavioral timing mode:
            [obj.b.params.xs, obj.b.params.xl, obj.b.params.xw] = deal(zeros(obj.b.params.n,obj.b.params.n));
            obj.defineUncertaintyKernelsWidths(false);
            for y = 1:obj.b.params.n 
                obj.b.params.xs(:,y) = normpdf(obj.b.params.t,y,obj.b.params.S(y))';           % small kernel, 
                obj.b.params.xl(:,y) = normpdf(obj.b.params.t,y,obj.b.params.L(y))';           % large kernel
                obj.b.params.xw(:,y) = normpdf(obj.b.params.t,y,obj.b.params.web(y))';         % Weber's law

            end
%             warning('I think we need to define L with the web somehow...the L depends on web...')
%             obj.b.params.xl = obj.b.params.xl.*obj.b.params.xw; % we will multiplex the weber frax with the large kernel
%             obj.b.params.xs = obj.b.params.xs;
            warning('changed this so we keep kernels going forward forever from the cue')
%             obj.b.params.xs(:,obj.b.params.oT)=0; 
%             obj.b.params.xl(:,obj.b.params.oT)=0; 
%             obj.b.params.xw(:,obj.b.params.oT)=0;     % leave out times outside trial
            obj.b.params.xs(:,1:obj.b.params.CS)=0; 
            obj.b.params.xl(:,1:obj.b.params.CS)=0; 
            obj.b.params.xw(:,1:obj.b.params.CS)=0;     % leave out before the cue

            obj.b.params.xs=obj.b.params.xs./sum(obj.b.params.xs); 
            obj.b.params.xl=obj.b.params.xl./sum(obj.b.params.xl); 
            obj.b.params.xw=obj.b.params.xw./sum(obj.b.params.xw);    % make prob dist's

            obj.b.params.xs(isnan(obj.b.params.xs))=0; 
            obj.b.params.xl(isnan(obj.b.params.xl))=0; 
            obj.b.params.xw(isnan(obj.b.params.xw))=0; % nan's to zeros

            if Plot
                [f1,ax1] = makeStandardFigure(3,[3,1]);
                title(ax1(1),'p(t|tau, \sigma_s) -- s means small width')
                title(ax1(2),'p(t|tau, \sigma_l * \sigma_{weber})')
                title(ax1(3),'p(t|tau, \sigma_{weber})')
                xlabel(ax1(3),'Veridical Time (s)')
                plot(ax1(1),obj.b.params.t_subjective,obj.b.params.xs)
                plot(ax1(2),obj.b.params.t_subjective,obj.b.params.xl)
                plot(ax1(3),obj.b.params.t_subjective,obj.b.params.xw)
                set(f1, 'userdata', ['defineUncertaintyKernels(obj,Plot) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
                set(f1, 'name', 'Behavioral timing uncertainty kernels')
            end 
        end
        function defineBeta(obj, Plot)
            if nargin < 2
                Plot=false;
            end
            % correction term:  [without correction | with correction]
            % this term updates the estimate of the value of veridical state t (Vh_t)
            % based on the current estimate of Vh_t + learning rate applied to the RPE
            % at subjective state tau (times probability of being in veridical t given
            % tau) - corrected (beta) estimate of value at tau times prob of being in t
            % given tau
            obj.b.params.beta = [zeros(obj.b.params.n,1), obj.params.alpha*(exp((log(obj.params.gamma))^2*(obj.b.params.L.^2-obj.b.params.S.^2)'/2)-1)];
            % the correction will be applied anytime sensory feedback is provided

            % we assume that the animal is constantly resolving its own state
            % uncertainty by virtue of the fact it must act upon its estimate of time
            if Plot
                [f1,ax1] = makeStandardFigure(2,[2,1]);
                title(ax1(1),'beta - state uncertainty, but no feedback')
                title(ax1(2),'beta - state uncertainty, + feedback')
                xlabel(ax1(2),'Veridical Time (s)')
                plot(ax1(1),obj.b.params.t_subjective,obj.b.params.beta(:,1))
                plot(ax1(2),obj.b.params.t_subjective,obj.b.params.beta(:,2))
                set(f1, 'userdata', ['defineBeta(obj, Plot) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
                set(f1, 'name', 'Behavioral timing betas')
            end 
        end
        function learning(obj, Mode, Plot, numIter)
            if nargin < 2
                Mode = 'uncertain+feedback'; % 'uncertain, no feedback', 'weber'
            end
            if nargin < 3
                Plot=true;
            end
            if nargin < 4
                numIter=1000;
            end
            % behavioral first...
            V = obj.b.params.V;
            T = obj.b.params.T;
            xl = obj.b.params.xl;
            xs = obj.b.params.xs;
            xw = obj.b.params.xw;
            r = obj.b.params.r;
            beta = obj.b.params.beta;
            gamma = obj.params.gamma;
            alpha = obj.params.alpha;
            n = obj.b.params.n;
            t_subjective = obj.b.params.t_subjective;
            
            % learning with feedback
            [w,... V_tau, the estimated value at subjective time tau
                Vh,... the estimate of the value function at veridical time t
                delta,... the rpe at subjective timing state tau
                ] = deal(zeros(n,3));              
            %
            % NB! w = Vh_tau, the subjective estimate of value. it's NOT the weight of
            % the state
            % 
            if Plot
                [f1,ax1] = makeStandardFigure(2,[2,1]);
                C = linspecer(10); 
                set(ax1(1), 'ColorOrder',C);%[blacks;C]);
                set(ax1(2), 'ColorOrder',C);
                title(ax1(2),'RPE')
                xlabel(ax1(2),'Veridical Time')
                legend(ax1(1),'show','Location','Northwest')
                set(f1, 'userdata', ['learning(obj, Mode, Plot, numIter) Mode=' Mode, ' numIter=', num2str(numIter), ' CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            end
            
            if strcmp(Mode, 'uncertain, no feedback') || Mode==1
                title(ax1(1),['$\hat{V_t}$, beta=', num2str(beta(1,1))], 'interpreter', 'latex')
                set(f1, 'name', 'State uncertainty without feedback, beta=0')
                %beta=0, no feedback
                warning('changed this search: ********************')
                for iter = 1:numIter % for each iteration of learning...
                    for y = 1:T+5%T+1 ********************
                        Vh(y,1) = w(:,1)'*xs(:,y); % Vh=Vt, w=Vtau, xs=p(t|tau, sigma=s) -- assuming we got feedback
                        Vh(y+1,1) = w(:,1)'*xl(:,y+1); % Vh=Vt+1, w=Vtau+1, xs=p(t+1|tau+1, sigma=l)
                        delta(y,1) = r(y) + gamma*Vh(y+1,1) - Vh(y,1);
                        w(:,1) = w(:,1) + (alpha*delta(y,1)-beta(y,1)*w(:,1)).*xs(:,y); 
%                         w(T+1:end,1) = 1;          % value stays high
%                         until reward ********************
                    end
                    if Plot && ismember(iter,(numIter/10).*(1:10))
                        Vhplot = Vh(:,1)./max(Vh(1:T,1)); % ********************
                        plot(ax1(1),t_subjective,Vhplot, 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
%                         plot(ax1(1),t_subjective,Vh(:,1), 'linewidth', 3,
%                         'displayname', ['iter=', num2str(iter)]) ********************
                        plot(ax1(2),t_subjective,delta(:,1), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
                    end
                    obj.b.unof.numIter = numIter;
                    obj.b.unof.w = w;
                    obj.b.unof.Vh = Vh;
                    obj.b.unof.delta = delta;
                end
            elseif strcmp(Mode, 'weber') || Mode==3
                set(f1, 'name', 'Weber uncertainty without feedback')
                title(ax1(1),'$\hat{V_t}$, no feedback (weber uncertainty)', 'interpreter', 'latex')
                warning('changed this search: ********************')
                for iter = 1:numIter
                    for y = 1:T+5%T+1 ********************
                        Vh(y,3) = w(:,3)'*xw(:,y);
                        Vh(y+1,3) = w(:,3)'*xw(:,y+1); % actually is more uncertain in the next state
                        delta(y,3) = r(y) + gamma*Vh(y+1,3) - Vh(y,3);
                        w(:,3) = w(:,3) + alpha*delta(y,3).*xw(:,y);
%                         w(T+1:end,3) = r(T);         	% value stays high
%                         until reward ********************
                    end
                    if ismember(iter,(numIter/10).*(1:10))
                        Vhplot = Vh(:,3)./max(Vh(1:T,3)); % ********************
                        plot(ax1(1),t_subjective,Vhplot, 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
%                         plot(ax1(1),t_subjective,Vh(:,3), 'linewidth', 3,
%                         'displayname', ['iter=', num2str(iter)]) ********************
                        plot(ax1(2),t_subjective,delta(:,3), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
                    end
                end
                obj.b.weber.numIter = numIter;
                obj.b.weber.w = w;
                obj.b.weber.Vh = Vh;
                obj.b.weber.delta = delta;
            else
                % beta=beta, with feedback
                title(ax1(1),['$\hat{V_t}$, beta=', num2str(beta(1,2))], 'interpreter', 'latex')
                set(f1, 'name', ['State uncertainty with feedback, beta=', num2str(beta(1,2))])
                warning('changed this search: ********************')
                for iter = 1:numIter % for each iteration of learning...
                    for y = 1:T+5%T+1 ********************
                        Vh(y,2) = w(:,2)'*xs(:,y);
                        Vh(y+1,2) = w(:,2)'*xl(:,y+1);
                        delta(y,2) = r(y) + gamma*Vh(y+1,2) - Vh(y,2);
                        w(:,2) = w(:,2) + (alpha*delta(y,2)-beta(y,2)*w(:,2)).*xs(:,y);
%                         w(T+1:end,2) = 1;          % value stays high
%                         until reward ********************
                    end
                    if ismember(iter,(numIter/10).*(1:10))
                        Vhplot = Vh(:,2)./max(Vh(1:T,2)); % ********************
                        plot(ax1(1),t_subjective,Vhplot, 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
%                         plot(ax1(1),t_subjective,Vh(:,2), 'linewidth', 3,
%                         'displayname', ['iter=', num2str(iter)]) ********************
                        plot(ax1(2),t_subjective,delta(:,2), 'linewidth', 3, 'displayname', ['iter=', num2str(iter)])
                    end
                end
                obj.b.uwithf.numIter = numIter;
                obj.b.uwithf.w = w;
                obj.b.uwithf.Vh = Vh;
                obj.b.uwithf.delta = delta;
            end
            plot(ax1(1),t_subjective,V, 'k', 'linewidth', 3, 'displayname', 'true V_t')
        end
        function plotJMFig1(obj)
            warning('must use learning(obj, Mode, Plot, numIter) in Mode 1 and 2 first')
            V = obj.b.params.V;
            T = obj.b.params.T;
            xl = obj.b.params.xl;
            xs = obj.b.params.xs;
            xw = obj.b.params.xw;
            r = obj.b.params.r;
            beta = obj.b.params.beta;
            gamma = obj.params.gamma;
            state_width_ms = obj.b.params.state_width_ms;
            n = obj.b.params.n;
            t_subjective = obj.b.params.t_subjective;
            Vh = [obj.b.unof.Vh(:,1),obj.b.uwithf.Vh(:,2),obj.b.weber.Vh(:,3)];
            
            [f1,ax2] = makeStandardFigure(2,[1,2]);
            set(f1, 'name', ['Weber and Feedback Models, beta=', num2str(beta(1,2))])
            set(f1, 'userdata', ['plotJMFig1(obj) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            fL = [3 1];                             % index of relevant value estimates
            hL = round([0.25*3300/state_width_ms 0.5*3300/state_width_ms 0.75*3300/state_width_ms]);                        % initial location of red curves
            y = round(0.25*3300/state_width_ms);                                 % length of red curves along x-axis
            cols=[205 92 92; 255 0 0; 150 0 0]/255; % red color schemes
            feed = {'Weber uncertainty, no feedback',['State uncertainty with continual feedback (beta=', num2str(beta(1,2)) ')']};
            for e = [2 1]
                f = fL(e);
                subplot(1,2,e) 
                plot(t_subjective,V, 'linewidth', 3)
                plot(t_subjective,Vh(:,f),'k', 'linewidth', 3)
                for g = 1:3
                    redcurve = nan(1,n);
            %         redcurve = nan(1,y);
                    h = hL(g);
                    redcurve(h+(0:y)) = Vh(h,f).*gamma.^(0:-1:-y); % thi is the estimate before feedback
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
        end
        function simulateBehavior(obj, varargin)
            pp = inputParser;
            addParameter(pp, 'Plot', true, @islogical);
            addParameter(pp, 'p_veridical', [0,1], @isnumeric);
            addParameter(pp, 'p_slow', [0.008,0.008], @isnumeric); %[0.0001,0.05], @isnumeric);
            addParameter(pp, 'p_fast', [0.0275,0.0275], @isnumeric);%[0.004,0.90], @isnumeric);
            addParameter(pp, 'n_sim', 100, @isnumeric);%[0.004,0.90], @isnumeric);
            parse(pp, varargin{:});
            Plot            = pp.Results.Plot;
            p_veridical 	= pp.Results.p_veridical;
            p_slow 			= pp.Results.p_slow;
            p_fast          = pp.Results.p_fast;
            n_sim           = pp.Results.n_sim;
            
            % simulate veridical timing behavior... state transitions
            % happen every obj.b.params.state_width_ms
            n = obj.b.params.n;
            obj.b.params.n_sim = n_sim;
            
            % simulate many times to get a CI
            tvos = nan(n_sim,n);
            tvss = nan(n_sim,n);
            tvfs = nan(n_sim,n);
            for i = 1:1000
                 tvos(i,:) = obj.runBehaviorSim(p_veridical);
                 tvss(i,:) = obj.runBehaviorSim(p_slow);
                 tvfs(i,:) = obj.runBehaviorSim(p_fast);
            end
            t_veridical_optimal = mean(tvos,1);
            t_veridical_slow = mean(tvss,1);
            t_veridical_fast = mean(tvfs,1);
            
            alph = 0.05;
            CI_l_tvos = sort(tvos,1); CI_l_tvos = CI_l_tvos(round(size(CI_l_tvos,1)*alph/2), :);
            CI_u_tvos = sort(tvos,1); CI_u_tvos= CI_u_tvos(round(size(CI_u_tvos,1)*(1-alph/2)), :);
            CI_l_tvss = sort(tvss,1); CI_l_tvss = CI_l_tvss(round(size(CI_l_tvss,1)*alph/2), :);
            CI_u_tvss = sort(tvss,1); CI_u_tvss = CI_u_tvss(round(size(CI_u_tvss,1)*(1-alph/2)), :);
            CI_l_tvfs = sort(tvfs,1); CI_l_tvfs = CI_l_tvfs(round(size(CI_l_tvfs,1)*alph/2), :);
            CI_u_tvfs = sort(tvfs,1); CI_u_tvfs = CI_u_tvfs(round(size(CI_u_tvfs,1)*(1-alph/2)), :);
            
            obj.b.sim.p_veridical = p_veridical;
            obj.b.sim.p_slow = p_slow; 	
            obj.b.sim.p_fast = p_fast; 
            obj.b.sim.t_veridical_optimal = t_veridical_optimal;
            obj.b.sim.t_veridical_slow = t_veridical_slow;
            obj.b.sim.t_veridical_fast = t_veridical_fast;
            obj.b.sim.CI_l_tvos = CI_l_tvos;
            obj.b.sim.CI_u_tvos = CI_u_tvos;
            obj.b.sim.CI_l_tvss = CI_l_tvss;
            obj.b.sim.CI_u_tvss = CI_u_tvss;
            obj.b.sim.CI_l_tvfs = CI_l_tvfs;
            obj.b.sim.CI_u_tvfs = CI_u_tvfs;    
            
            if Plot
                obj.plotBehaviorSim;
            end
        end
        function plotBehaviorSim(obj)
            p_veridical = obj.b.sim.p_veridical;
            p_slow = obj.b.sim.p_slow; 	
            p_fast = obj.b.sim.p_fast; 
            t_veridical_optimal = obj.b.sim.t_veridical_optimal;
            t_veridical_slow = obj.b.sim.t_veridical_slow;
            t_veridical_fast = obj.b.sim.t_veridical_fast;
            CI_l_tvos = obj.b.sim.CI_l_tvos;
            CI_u_tvos = obj.b.sim.CI_u_tvos;
            CI_l_tvss = obj.b.sim.CI_l_tvss;
            CI_u_tvss = obj.b.sim.CI_u_tvss;
            CI_l_tvfs = obj.b.sim.CI_l_tvfs;
            CI_u_tvfs = obj.b.sim.CI_u_tvfs;  


            T = obj.b.params.T;
            beta = obj.b.params.beta;
            V = obj.b.params.V; V = V(1:T);
            Vh = obj.b.uwithf.Vh; Vh = Vh(1:T,2);
            delta = obj.b.uwithf.delta; delta = delta(1:T,2);
            Fs = 100;
            [b,a] = butter(2, 10 / (Fs/2));
            smoothdelta = filter(b,a,delta); 
            delta=delta(1:end-2);% trim off the edge artifact -- times after first-lick not defined
            smoothdelta = smoothdelta(1:end-2);

            
            t_subjective = obj.b.params.t_subjective;
            T = obj.b.params.T;
            
            
            [f1,ax1] = makeStandardFigure();
            [f2,ax2] = makeStandardFigure(3,[1,3]);
            title(ax1(1),'simulated subjective state space')
            ylabel(ax1(1),'Subjective Time')
            xlabel(ax1(1),'Veridical Time')
            

            set(f1, 'name', ['State Uncertainty at different subjective pacemaker rates, beta=', num2str(beta(1,2)), ' nsim=' num2str(obj.b.params.n_sim)]);
            set(f2, 'units', 'normalized', 'position', [0.1,0.2,0.6,0.25]);
            set(f2, 'name', ['State Uncertainty + Feedback Models at different subjective pacemaker rates, beta=', num2str(beta(1,2)), ' nsim=' num2str(obj.b.params.n_sim)]);

            set(f1, 'userdata', ['obj.simulateBehavior() CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            plot(ax1(1), [min(t_veridical_slow), max(t_veridical_slow)],[0,0],  'c--', 'linewidth', 2, 'displayname', 'cue')
            plot(ax1(1), [0,0],[min(t_subjective), max(t_subjective)],  'c--', 'linewidth', 2, 'displayname', 'cue', 'HandleVisibility', 'off')
            plot(ax1(1), [min(t_veridical_slow), max(t_veridical_slow)],[3.3,3.3],  'g--', 'linewidth', 2, 'displayname', 'first-lick')
            plot(ax1(1), [t_veridical_fast(T)-(t_subjective(T)-3.3),t_veridical_fast(T)-(t_subjective(T)-3.3)], [min(t_subjective), max(t_subjective)], 'r--', 'linewidth', 1, 'displayname', 'first-lick', 'HandleVisibility', 'off')
            plot(ax1(1), [t_veridical_optimal(T)-(t_subjective(T)-3.3),t_veridical_optimal(T)-(t_subjective(T)-3.3)], [min(t_subjective), max(t_subjective)], 'k--', 'linewidth', 1, 'displayname', 'first-lick', 'HandleVisibility', 'off')
            plot(ax1(1), [t_veridical_slow(T)-(t_subjective(T)-3.3),t_veridical_slow(T)-(t_subjective(T)-3.3)], [min(t_subjective), max(t_subjective)], 'b--', 'linewidth', 1, 'displayname', 'first-lick', 'HandleVisibility', 'off')
            
            plot(ax1(1), CI_l_tvfs,t_subjective,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax1(1), CI_u_tvfs,t_subjective,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax1(1), t_veridical_fast,t_subjective,  'b-', 'linewidth', 3, 'displayname', ['fast p[!stay, txn]=', mat2str(p_fast)])          
            
            plot(ax1(1), CI_l_tvos,t_subjective,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax1(1), CI_u_tvos,t_subjective,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax1(1), t_veridical_optimal,t_subjective,  'k-', 'linewidth', 3, 'displayname', ['optimal p[!stay, txn]=', mat2str(p_veridical)])
           
            plot(ax1(1), CI_l_tvss,t_subjective,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax1(1), CI_u_tvss,t_subjective,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax1(1), t_veridical_slow,t_subjective,  'r-', 'linewidth', 3, 'displayname', ['slow p[!stay, txn]=', mat2str(p_slow)])
            
            legend(ax1(1),'show','Location','best')



            set(f2, 'userdata', ['obj.simulateBehavior() CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])

         
            plot(ax2(1), [0,0],[0, max(Vh)],  'c--', 'linewidth', 2, 'displayname', 'cue')
            plot(ax2(1), [min(t_veridical_slow), max(CI_u_tvss)],[1,1],  'g--', 'linewidth', 2, 'displayname', 'threshold')
            plot(ax2(1), CI_l_tvfs(1:numel(Vh)),Vh,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(1), CI_u_tvfs(1:numel(Vh)),Vh,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(1), t_veridical_fast(1:numel(Vh)),Vh,  'b-', 'linewidth', 3, 'displayname', ['fast p[!stay, txn]=', mat2str(p_fast)])          
            plot(ax2(1), CI_l_tvos(1:numel(Vh)),Vh,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(1), CI_u_tvos(1:numel(Vh)),Vh,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(1), t_veridical_optimal(1:numel(Vh)),Vh,  'k-', 'linewidth', 3, 'displayname', ['optimal p[!stay, txn]=', mat2str(p_veridical)])
            plot(ax2(1), CI_l_tvss(1:numel(Vh)),Vh,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(1), CI_u_tvss(1:numel(Vh)),Vh,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(1), t_veridical_slow(1:numel(Vh)),Vh,  'r-', 'linewidth', 3, 'displayname', ['slow p[!stay, txn]=', mat2str(p_slow)])
            legend(ax2(1),'show','Location','best')
            
            plot(ax2(2), [0,0],[0, max(CI_u_tvss)],  'c--', 'linewidth', 2, 'displayname', 'cue')
            plot(ax2(2), [min(CI_u_tvss),max(CI_u_tvss)],[max(delta),max(delta)],  'g--', 'linewidth', 2, 'displayname', 'threshold')
            plot(ax2(2), CI_l_tvfs(1:numel(delta)),delta,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(2), CI_u_tvfs(1:numel(delta)),delta,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(2), t_veridical_fast(1:numel(delta)),delta,  'b-', 'linewidth', 3, 'displayname', ['fast p[!stay, txn]=', mat2str(p_fast)])          
            plot(ax2(2), CI_l_tvos(1:numel(delta)),delta,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(2), CI_u_tvos(1:numel(delta)),delta,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(2), t_veridical_optimal(1:numel(delta)),delta,  'k-', 'linewidth', 3, 'displayname', ['optimal p[!stay, txn]=', mat2str(p_veridical)])
            plot(ax2(2), CI_l_tvss(1:numel(delta)),delta,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(2), CI_u_tvss(1:numel(delta)),delta,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(2), t_veridical_slow(1:numel(delta)),delta,  'r-', 'linewidth', 3, 'displayname', ['slow p[!stay, txn]=', mat2str(p_slow)])
            ylim(ax2(2),[0,max(delta)])


            plot(ax2(3), [0,0], [0, max(CI_u_tvss)],  'c--', 'linewidth', 2, 'displayname', 'cue')
            plot(ax2(3), [0, max(CI_u_tvss)],[max(smoothdelta),max(smoothdelta)],  'g--', 'linewidth', 2, 'displayname', 'threshold')
            plot(ax2(3), CI_l_tvfs(1:numel(smoothdelta)),smoothdelta,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(3), CI_u_tvfs(1:numel(smoothdelta)),smoothdelta,  'b-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(3), t_veridical_fast(1:numel(smoothdelta)),smoothdelta,  'b-', 'linewidth', 3, 'displayname', ['fast p[!stay, txn]=', mat2str(p_fast)])          
            plot(ax2(3), CI_l_tvos(1:numel(smoothdelta)),smoothdelta,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(3), CI_u_tvos(1:numel(smoothdelta)),smoothdelta,  'k-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(3), t_veridical_optimal(1:numel(smoothdelta)),smoothdelta,  'k-', 'linewidth', 3, 'displayname', ['optimal p[!stay, txn]=', mat2str(p_veridical)])
            plot(ax2(3), CI_l_tvss(1:numel(smoothdelta)),smoothdelta,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')          
            plot(ax2(3), CI_u_tvss(1:numel(smoothdelta)),smoothdelta,  'r-', 'linewidth', 1, 'HandleVisibility', 'off')
            plot(ax2(3), t_veridical_slow(1:numel(smoothdelta)),smoothdelta,  'r-', 'linewidth', 3, 'displayname', ['slow p[!stay, txn]=', mat2str(p_slow)])
            ylim(ax2(3),[0,max(smoothdelta)])
                        
            ylabel(ax2(1),'Subjective Value, V_τ = V_t ')
            ylabel(ax2(2),'dV_τ / dτ  ≈ RPE, δ_τ')
            ylabel(ax2(3),'Smooth RPE, δ_τ')

            for ii = 1:3
                xlim(ax2(ii), [-0.5,5.5])
                xticks(ax2(ii), -0.5:0.5:7.5)
                xlabel(ax2(ii), 'Veridical Time, t')
            end
            
        end
        function t_veridical = runBehaviorSim(obj, p)
            % p = [1-p_stay (should be the case for times <
            % obj.b.params.state_width_ms for veridical timer), 
            % p_txn (should be 1 at obj.b.params.state_width_ms)]
            n = obj.b.params.n;
            CS = obj.b.params.CS;
            t_subjective = obj.b.params.t_subjective;
            T_time = obj.b.params.T_time;
            
            % don't bother modeling pre-cue
            t_veridical = t_subjective(t_subjective<=0);
            t_ms = 0; % relative to cue
            
            t_breaks = obj.b.params.state_width_ms;
            t_since_break = 0;
            
            for t = 0:1000*T_time+2000%1000
                t_ms = t_ms+1;
                t_since_break = t_since_break+1;
                if numel(t_veridical) == numel(t_subjective)
                    break
                end
                if t_since_break == t_breaks
                    % should transition with high probability
                    while true
                        if obj.flip(p(2))
                            t_veridical(end+1) = t_ms/1000;
                            t_since_break = 0;
                            break
                        else
                            t_ms = t_ms+1;
                        end
                    end
                else
                    % should transition with low probability
                    if obj.flip(p(1))
                        t_veridical(end+1) = t_ms/1000;
                        continue
                    end
                end
            end
        end
        function r = flip(obj,p)
            %
            % if p is prob, is a weighted coin tos
            % if p is vector, randomly picks one from the vector
            %
            if numel(p) == 1
                r = rand < p;
            else
                r = p(randi(numel(p)));
            end
        end
        function simBehaviorValueTrajectory(obj)
            %
            %   NB: plotsimulation is preferred, as it shows CIs
            %

            % must call obj.simulateBehavior first. also must have learned
            % value trajectory already learning(obj, Mode, Plot, numIter) in Mode 2
            T = obj.b.params.T;
            beta = obj.b.params.beta;
            V = obj.b.params.V; V = V(1:T);
            Vh = obj.b.uwithf.Vh; Vh = Vh(1:T,2);
            delta = obj.b.uwithf.delta; delta = delta(1:T,2);

            t_subjective = obj.b.params.t_subjective; t_subjective = t_subjective(1:T);
            
            t_veridical_fast = obj.b.sim.t_veridical_fast; t_veridical_fast = t_veridical_fast(1:T);
            t_veridical_optimal = obj.b.sim.t_veridical_optimal; t_veridical_optimal = t_veridical_optimal(1:T);
            t_veridical_slow = obj.b.sim.t_veridical_slow; t_veridical_slow = t_veridical_slow(1:T);

            % filter params for simulating perceptual lag of cue and calcium indicator dynamics
            Fs = 100;
            [b,a] = butter(2, 10 / (Fs/2));
            smoothdelta = filter(b,a,delta);

            [f1,ax1] = makeStandardFigure(3,[1,3]);
            set(f1, 'units', 'normalized', 'position', [0.1,0.2,0.6,0.25])
            set(f1, 'name', ['Feedback Models at different subjective pacemaker rates, beta=', num2str(beta(1,2))])
            set(f1, 'userdata', ['simBehaviorValueTrajectory(obj) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            title(ax1(1), 'Estimated Value')
            title(ax1(2), 'RPE')
            title(ax1(3), 'Smoothed RPE')
            
            plot(ax1(1),t_subjective,V,'g-','linewidth', 4, 'displayname', 'True Value')
            plot(ax1(1),t_veridical_fast,Vh,'b-','linewidth', 3, 'displayname', 'Vh, fast pacemaker')
            plot(ax1(1),t_veridical_optimal,Vh,'k-','linewidth', 3, 'displayname', 'Vh, optimal timer')
            plot(ax1(1),t_veridical_slow,Vh,'r-','linewidth', 3, 'displayname', 'Vh, slow pacemaker')
            
            
            plot(ax1(2),t_veridical_fast(1:end-2),delta(1:end-2),'b-','linewidth', 3, 'displayname', 'Vh, fast pacemaker') % truncated to avoid edge artifact
            plot(ax1(2),t_veridical_optimal(1:end-2),delta(1:end-2),'k-','linewidth', 3, 'displayname', 'Vh, optimal timer')
            plot(ax1(2),t_veridical_slow(1:end-2),delta(1:end-2),'r-','linewidth', 3, 'displayname', 'Vh, slow pacemaker')

            plot(ax1(3),t_veridical_fast(1:end-2),smoothdelta(1:end-2),'b-','linewidth', 3, 'displayname', 'Vh, fast pacemaker') % truncated to avoid edge artifact
            plot(ax1(3),t_veridical_optimal(1:end-2),smoothdelta(1:end-2),'k-','linewidth', 3, 'displayname', 'Vh, optimal timer')
            plot(ax1(3),t_veridical_slow(1:end-2),smoothdelta(1:end-2),'r-','linewidth', 3, 'displayname', 'Vh, slow pacemaker')
            
            ylabel(ax1(1),'Value')
            ylim(ax1(1),[0 1])
            xlim(ax1(1),[min(t_veridical_slow) max(t_veridical_slow)])
            
            ylabel(ax1(2),'RPE')
            xlim(ax1(2),[min(t_veridical_slow) max(t_veridical_slow)])

            ylabel(ax1(3),'smoothed RPE')
            xlim(ax1(3),[min(t_veridical_slow) max(t_veridical_slow)])

            
            legend('show','location', 'best', 'box','off')
        end
        
        
        
        function learnPerceptualTask(obj, varargin)
            pp = inputParser;
            addParameter(pp, 'Plot', true, @islogical);
            addParameter(pp, 'p_veridical', [0,1], @isnumeric);
            addParameter(pp, 'p_slow', [0.0001,0.05], @isnumeric);
            addParameter(pp, 'p_fast', [0.004,0.90], @isnumeric);
            addParameter(pp, 'ntrials', 1000, @isnumeric);
            addParameter(pp, 'nofeedback', false, @islogical);
            addParameter(pp, 'prob_correct_mode', false, @islogical);
            parse(pp, varargin{:});
            Plot            = pp.Results.Plot;
            p_veridical 	= pp.Results.p_veridical;
            p_slow 			= pp.Results.p_slow;
            p_fast          = pp.Results.p_fast;
            ntrials         = pp.Results.ntrials;
            nofeedback      = pp.Results.nofeedback;
            prob_correct_mode = pp.Results.prob_correct_mode;
            
            if nofeedback
                warning('nofeedback')
                obj.p.nofeedback = true;
            else
                obj.p.nofeedback = false;
            end
            if prob_correct_mode
                warning('prob_correct_mode')
                obj.p.prob_correct_mode = true;
            else
                obj.p.prob_correct_mode = false;
            end
            
            
            [f1,ax1] = makeStandardFigure(3,[3,1]);
            C = linspecer(10); 
            set(ax1(1), 'ColorOrder',C);%[blacks;C]);
            set(ax1(2), 'ColorOrder',C);
            set(ax1(3), 'ColorOrder',C);
            title(ax1(2),'RPE')
            title(ax1(3),'Vtau (w)')
            xlabel(ax1(2),'Veridical Time')
            legend(ax1(1),'show','Location','Northwest')
            set(f1, 'userdata', ['obj.learnPerceptualTask() CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            
            % simulate perceptual timing behavior... state transitions
            % happen every obj.b.params.state_width_ms
            gamma = obj.params.gamma;
            alpha = obj.params.alpha;
            
%             for trial
%                 draw a time
%                 progress through time in subjective space
%                 backpropagate the value

            % for the purposes of learning, we will assume the animal is an
            % optimal timer, since we are uncertain of the time of the
            % trial that's probs enough to deal with for now
            
            state_width_ms = 100;
%             n = round((3000-1)/state_width_ms);
            n = round((3000+400)/state_width_ms)+1;
            t_subjective = linspace(-0.4,3,n);
            w = 0.15;
            web = w*(1:n);
            S = .1+zeros(1,n);
            oT = find(t_subjective <= 0);
            
%             xw = zeros(n,n);xs = zeros(n,n);
%             for y = 1:n
%                 xw(:,y) = normpdf(1:n,y,web(y))';
%                 xs(:,y) = normpdf(1:n,y,S(y))';
%             end
%             xw(:,oT)=0; xs(:,oT)=0;     % leave out times outside trial
%             xw=xw./sum(xw); xs=xs./sum(xs);    % make prob dist's
%             xw(isnan(xw))=0; xs(isnan(xs))=0; % nan's to zeros
            xw = [];
            xs = [];
            
            
            
            warning('Debug:reward delivered/omitted at second tone')
            Ts = [find(t_subjective>=0.6,1,'first'), find(t_subjective>=1.05,1,'first'), find(t_subjective>=1.26,1,'first'), find(t_subjective>=1.74,1,'first'), find(t_subjective>=1.95,1,'first'), find(t_subjective>=2.4,1,'first')];
            n_trials = zeros(6,1);
            n_complete = zeros(6,1);
            n_rewards = zeros(6,1);
            
            pbail = 1/(2*n);
%             pbail = 0;
            warning('Debug: flipping only a few')
%             ntrials = 5000;
            [w,... V_tau, the estimated value at subjective time tau
                Vh,... the estimate of the value function at veridical time t
                delta,... the rpe at subjective timing state tau
                ] = deal(zeros(n,1));
            for trial = 1:ntrials
%                 if trial==1
%                     ix = 6;
%                 else
%                     ix = obj.flip(1:6);
                      ix=obj.flip([1,2,3,4,5,6]);
%                 end
                T = Ts(ix);
%                 T = Ts(4);
                % assign reward based on how likely to guess correctly at
                % each for now
                r = zeros(n,1); 
                if ix==1, if obj.flip(0.98),r(T) = 1; else,r(T) = 0;end
                elseif ix==2, if obj.flip(0.94), r(T) = 1; else, r(T) = 0;end
                elseif ix==3, if obj.flip(0.82), r(T) = 1; else, r(T) = 0;end
                elseif ix==4, if obj.flip(0.71), r(T) = 1; else, r(T) = 0;end
                elseif ix==5, if obj.flip(0.85), r(T) = 1; else, r(T) = 0;end
                elseif ix==6, if obj.flip(0.94), r(T) = 1; else, r(T) = 0;end
                end
                n_trials(ix) = n_trials(ix)+1;
                [Vh,delta,w,n_complete,n_rewards,rewarded] = obj.tryWeberUpdate(T, Vh, delta, w, trial,n_complete,n_rewards,pbail,xw,xs,r,gamma,alpha,ix,t_subjective,obj.p.nofeedback,obj.p.prob_correct_mode);
                if ismember(trial,(ntrials/10).*[1:10])
                    obj.plotLearningPerception(ax1, t_subjective, Vh, delta, w, trial, rewarded);
                end
            end
            plot(ax1(1),t_subjective,Vh, 'k-', 'linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
            plot(ax1(2),t_subjective,delta, 'k-','linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
            plot(ax1(3),t_subjective,w, 'k-','linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
            if trial > 1000
                legend(ax1(1),'off')
            end
            
            [f1,ax1] = makeStandardFigure(1,[1,1]);
            bar(ax1, n_trials, 'k', 'displayname', 'number of trials initialized in each category')
            bar(ax1, n_complete, 'displayname', 'number of trials completed in each category')
            bar(ax1, n_rewards, 'g', 'displayname', 'number of trials rewarded in each category')
            legend('show')
            
            obj.p.weber.ntrials = ntrials;
            obj.p.weber.n = n;
            obj.p.weber.w = w;
            obj.p.weber.Vh = Vh;
            obj.p.weber.delta = delta;
            obj.p.weber.n_trials = n_trials;
            obj.p.weber.n_complete = n_complete;
            obj.p.weber.n_rewards = n_rewards;
            obj.p.weber.t_subjective = t_subjective;
            obj.p.weber.Ts = Ts;
            obj.p.weber.pbail = pbail;
            obj.p.weber.xw = xw;
            obj.p.weber.xs = xs;
            
        end
        function plotLearningPerception(obj, ax1, t_subjective, Vh, delta, w, trial, rewarded)
            if rewarded(1) == -1
                title(ax1(1), ['Last trial bailed @ T=', num2str(rewarded(2)) 's'])
            elseif rewarded(1) == 1
                title(ax1(1), ['Last trial was rewarded, T=', num2str(rewarded(2))])
            elseif rewarded(1) == 0
                title(ax1(1), ['Last trial was NOT rewarded, T=', num2str(rewarded(2))])
            else
                title(ax1(1), '??')
            end
            plot(ax1(1),t_subjective,Vh, 'linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
            plot(ax1(2),t_subjective,delta, 'linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
            xlim(ax1(2), get(ax1(1),'xlim'));
            plot(ax1(3),t_subjective,w, 'linewidth', 3, 'displayname', ['trial#=', num2str(trial)])
        end
        function [Vh,delta,w,n_complete,n_rewards, rewarded] = tryWeberUpdate(obj,T, Vh, delta, w, trial,n_complete,n_rewards,pbail,xw,xs,r,gamma,alpha,ix,t_subjective,nofeedback,prob_correct_mode)
            % if the animal bails, exit without updating...
            
            
            rewarded = [-1,t_subjective(T)];
            Vh_init=Vh; delta_init=delta; w_init=w;
            %
            % Debug: let's also assume a negative RPE if we bail out on the
            % task...
            %
            bailnegativeRPE = true;
            thistrialbailed=false;
            if bailnegativeRPE
                for t=1:T
                    if obj.flip(pbail)
                        T=t;
                        r(T)=0;
                        thistrialbailed=true;
                        rewarded = [-1,t_subjective(T)];
                        break
                    end
                end
            end
                

            oT = find(t_subjective <= 0);
            n=numel(t_subjective);
            
            S = zeros(1,n);
            S(max(oT)+1:end) = 0.15.*(1:numel(Vh)-max(oT)); %0.15.*(1:numel(Vh)); 
            S(t_subjective <= 0) = 0;
            S(T+1) = 2*0.15;
            S(T+2:end) = 0;
%             S=S./sum(S);
            
            L = S;
            if ~nofeedback
                S(T) = 0.15;
                S(T+1) = 2*0.15;
            else
                S(T) = 0.15;
                L = S;
            end
            
            
            
            ss = zeros(n,n);ll = zeros(n,n);
            for y = 1:n
                ss(:,y) = normpdf(1:n,y,S(y))';
                ll(:,y) = normpdf(1:n,y,L(y))';
            end
            ss(:,oT)=0; ll(:,oT)=0;     % leave out times outside trial
            ss=ss./sum(ss); ll=ll./sum(ll);    % make prob dist's
            ss(isnan(ss))=0; ll(isnan(ll))=0; % nan's to zeros
            
            beta = [zeros(numel(Vh),1), alpha*(exp((log(gamma))^2*(L.^2-S.^2)'/2)-1)];
            
            delta = nan(1,n);
            %
            % Debug: let's assume when the animal gets 2nd tone, it does
            % the t+1 comparison with the feedback state... and let's
            % assume the value of the feedback state is proportional to
            % prob of getting correct, i.e., we will set the value
            % based on reward no reward there...
            %
            % assume mouse learns in next state if correct or
            % not... However, we won't carry over the reward/omission in
            % the learned value function because we are focused on times
            % before feedback in the value landscape
            p_correct = [0.98,0.94,0.82,0.71,0.85,0.94];
            if r(T)==1
%                 Vh(T+2:end) = 1;
                if prob_correct_mode
                    w(T+1:end) = 1;         	% value is 1 after reward
                else
                    w(T+1:end) = 1;         	% value is 1 after reward
                end
%                 disp('rewarded')
            else
%                 Vh(T+2:end) = 0;
                w(T+1:end) = 0;         	% value is 0 after mistake
%                 disp('NOT rewarded')
            end
            if prob_correct_mode && r(T)==1
                r(T) = p_correct(ix); 
            end
            
            
            
            
            for t = 1:T+1% T+1
                if ~bailnegativeRPE && obj.flip(pbail) && t < T
                    Vh=Vh_init;
                    w=w_init;
                    disp('bailed')
                    return
                end
                if t==T
                    if thistrialbailed
                    else
                        n_complete(ix) = n_complete(ix)+1;
                        if r(T)>0
                            n_rewards(ix) = n_rewards(ix)+1;
                            rewarded(1)=1;
                        else
                            rewarded(1)=0;
                        end
                    end
                end
                Vh(t) = w'*ss(:,t);
                Vh(t+1) = w'*ll(:,t+1); % actually is more uncertain in the next state
                delta(t) = r(t) + gamma*Vh(t+1) - Vh(t);
                if nofeedback
                    w = w + alpha*delta(t).*ss(:,t);
                else
                    w = w + (alpha*delta(t)-beta(t)*w).*ss(:,t);
                end
            end
            % before kicking it back, we need to take out the stamps after
            % the T

%             Vh(T+1:end) = Vh_init(T+1:end);
%             w(T+1:end) = w_init(T+1:end);
            Vh(T:end) = Vh_init(T:end);
            w(T:end) = w_init(T:end);
            
        end
        function test_RPE(obj,ix)
            ntrials = obj.p.weber.ntrials;
            w = obj.p.weber.w;
            Vh = obj.p.weber.Vh;
            n = numel(Vh);
            delta = obj.p.weber.delta;
            n_trials = obj.p.weber.n_trials;
            n_complete = obj.p.weber.n_complete;
            n_rewards = obj.p.weber.n_rewards;
            t_subjective = obj.p.weber.t_subjective;
            alpha = obj.params.alpha;
            gamma = obj.params.gamma;
            Ts = obj.p.weber.Ts;
            T = Ts(ix);
            trial = ntrials+1;
            pbail = obj.p.weber.pbail;
            xw = obj.p.weber.xw;
            
            r = zeros(n,1); 
            if ix==1 
                    if obj.flip(0.98)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                elseif ix==2 
                    if obj.flip(0.94)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                elseif ix==3 
                    if obj.flip(0.82)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                 elseif ix==4
                    if obj.flip(0.71)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                 elseif ix==5
                    if obj.flip(0.85)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                 elseif ix==6
                    if obj.flip(0.94)
                        r(T) = 1;
                    else
                        r(T) = 0;
                    end
                end
            disp(['rewarded? =', num2str(r(Ts(ix)))])
            
            [f1,ax1] = makeStandardFigure(6,[3,2]);
            C = linspecer(10); 
            set(ax1(1), 'ColorOrder',C);%[blacks;C]);
            set(ax1(2), 'ColorOrder',C);
            set(ax1(3), 'ColorOrder',C);
            title(ax1(1),'Learned')
            title(ax1(2),['Ts=',num2str(ix)])
            ylabel(ax1(1),'Vh')
            ylabel(ax1(3),'RPE')
            ylabel(ax1(5),'Vtau (w)')
            xlabel(ax1(2),'Veridical Time')
            legend(ax1(1),'show','Location','Northwest')
            set(f1, 'userdata', ['obj.test_RPE() CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            obj.plotLearningPerception(ax1([1,3,5]), t_subjective, Vh, delta, w, trial);
            
            [Vh,delta,w,~,~] = obj.tryWeberUpdate(Ts(ix), Vh, delta, w, trial,n_complete,n_rewards,pbail,xw,r,gamma,alpha,ix);
            obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh, delta, w, trial);
        end
        function test_RPE2(obj,ix,nextisrewarded,fast_slowMode)
            if nargin < 3
                nextisrewarded = true;
            end
            if nargin < 4
                fast_slowMode = true;
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
                    [Vh,delta,w,~,~,rewarded] = tryWeberUpdate(obj,T, Vh, delta, w, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,obj.p.nofeedback,obj.p.prob_correct_mode);
                    obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh, delta, w, trial,rewarded);
                elseif fast_slowMode
                    % I think we can model fast and slow clock by stretching Vh
                    % and V in the weber update
                    % I think this makes sense because the Vh represents the
                    % subjective value estimate?
                    %
                    if ix <= 3 % we need relatively slow clock to pick short time
                        [Vh_slow,delta,w_slow,~,~,rewarded] = tryWeberUpdate(obj,T, Vh_slow, delta, w_slow, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,obj.p.nofeedback,obj.p.prob_correct_mode);
                        obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh_slow, delta, w_slow, trial,rewarded);
                        
                    else % we need relatively fast clock to pick long time (assume the proper fxn so not making data up)
                        [Vh,delta,w,~,~,rewarded] = tryWeberUpdate(obj,T, Vh, delta, w, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,obj.p.nofeedback,obj.p.prob_correct_mode);
                        obj.plotLearningPerception(ax1([2,4,6]), t_subjective, Vh, delta, w, trial,rewarded);
                    end
                end
            end
            if fast_slowMode
                title(ax1(2),'short assumed to be slow clock')
            else
                title(ax1(2),['Ts=',num2str(ix)])
            end
        end
        function simulate_perception(obj,nsim)
            %
            % We will take the learned value landscape and simulate trials
            % of a mouse...
            %
            ntrials = obj.p.weber.ntrials;w = obj.p.weber.w;Vh = obj.p.weber.Vh;n = numel(Vh);delta = obj.p.weber.delta;n_complete = obj.p.weber.n_complete;n_rewards = obj.p.weber.n_rewards;t_subjective = obj.p.weber.t_subjective;alpha = obj.params.alpha;gamma = obj.params.gamma;Ts = obj.p.weber.Ts;
            pbail = obj.p.weber.pbail;
            
            [f1,ax1] = makeStandardFigure(3,[3,1]);
            C = linspecer(10); 
            set(ax1(1), 'ColorOrder',C);%[blacks;C]);
            set(ax1(2), 'ColorOrder',C);
            set(ax1(3), 'ColorOrder',C);
            title(ax1(2),'RPE')
            title(ax1(3),'Vtau (w)')
            xlabel(ax1(2),'Veridical Time')
            legend(ax1(1),'show','Location','Northwest')
            set(f1, 'userdata', ['obj.simulate_perception(nsim) CLASS_value_landscaping_obj runID=', num2str(obj.params.runID)])
            
            
            Value_correct = nan(6,n);
            Value_incorrect = nan(6,n);
            Value_bail = nan(6,n);
            RPE_correct = nan(6,n);
            RPE_incorrect = nan(6,n);
            RPE_bail = nan(6,n);
            
            Value_correct_all_1 = [];
            Value_incorrect_all_1 = [];
            RPE_correct_all_1 = [];
            RPE_incorrect_all_1 = [];
            
            Value_correct_all_2 = [];
            Value_incorrect_all_2 = [];
            RPE_correct_all_2 = [];
            RPE_incorrect_all_2 = [];
            
            Value_correct_all_3 = [];
            Value_incorrect_all_3 = [];
            RPE_correct_all_3 = [];
            RPE_incorrect_all_3 = [];
            
            Value_correct_all_4 = [];
            Value_incorrect_all_4 = [];
            RPE_correct_all_4 = [];
            RPE_incorrect_all_4 = [];
            
            Value_correct_all_5 = [];
            Value_incorrect_all_5 = [];
            RPE_correct_all_5 = [];
            RPE_incorrect_all_5 = [];
            
            Value_correct_all_6 = [];
            Value_incorrect_all_6 = [];
            RPE_correct_all_6 = [];
            RPE_incorrect_all_6 = [];
            
            Value_bail_all = nan(6,n);
            RPE_bail_all = nan(6,n);
            n_sims = zeros(1,6);
            n_bail = zeros(1,6);
            n_correct = zeros(1,6);
            n_incorrect = zeros(1,6);
            for trial = 1:nsim
                ix=obj.flip([1,2,3,4,5,6]);
                T = Ts(ix);
                r = zeros(n,1); 
                if ix==1, if obj.flip(0.98),r(T) = 1; else,r(T) = 0;end
                elseif ix==2, if obj.flip(0.94), r(T) = 1; else, r(T) = 0;end
                elseif ix==3, if obj.flip(0.82), r(T) = 1; else, r(T) = 0;end
                elseif ix==4, if obj.flip(0.71), r(T) = 1; else, r(T) = 0;end
                elseif ix==5, if obj.flip(0.85), r(T) = 1; else, r(T) = 0;end
                elseif ix==6, if obj.flip(0.94), r(T) = 1; else, r(T) = 0;end
                end
                n_sims(ix) = n_sims(ix)+1;
                [Vh,delta,w,n_complete,n_rewards,rewarded] = obj.tryWeberUpdate(T, Vh, delta, w, trial,n_complete,n_rewards,pbail,[],[],r,gamma,alpha,ix,t_subjective,obj.p.nofeedback,obj.p.prob_correct_mode);
                if rewarded(1) == 1
                    Value_correct(ix,:) = nanmean([Value_correct(ix,:); Vh']);
                    RPE_correct(ix,:) = nanmean([RPE_correct(ix,:); delta]);
                    if ix==1
                        Value_correct_all_1(end+1,:) = Value_correct(ix,:);
                        RPE_correct_all_1 = RPE_correct(ix,:);
                    end
                    n_correct(ix) = n_correct(ix) + 1;
                elseif rewarded(1) == 0
                    n_incorrect(ix) = n_incorrect(ix) + 1;
                    Value_incorrect(ix,:) = nanmean([Value_incorrect(ix,:); Vh']);
                    RPE_incorrect(ix,:) = nanmean([RPE_incorrect(ix,:); delta]);
                else
                    % bail
                    n_bail(ix) = n_bail(ix) + 1;
                    Value_bail(ix,:) = nanmean([Value_bail(ix,:); Vh']);
                    RPE_bail(ix,:) = nanmean([RPE_bail(ix,:); delta]);
                end
                if ismember(trial,(nsim/10).*[1:10])
                    obj.plotLearningPerception(ax1, t_subjective, Vh, delta, w, trial, rewarded);
                end
            end
            [f,ax]=plotperceptionsim(obj,[],[],Value_correct,RPE_correct,Value_incorrect,RPE_incorrect,Value_bail,RPE_bail,t_subjective, n_correct, n_incorrect, n_bail);
            
        end
        function [f,ax]=plotperceptionsim(obj,f,ax,Value_correct,RPE_correct,Value_incorrect,RPE_incorrect,Value_bail,RPE_bail,t_subjective, n_correct, n_incorrect, n_bail)
            if isempty(f)
                [f,ax]=makeStandardFigure(6,[2,3]);
            end
            title(ax(1),'Correct')
            title(ax(2),'Incorrect')
            title(ax(3),'Bail')
            ylabel(ax(1),'Mean Vh')
            ylabel(ax(4),'Mean \delta')

            plot(ax(1),t_subjective,Value_correct')
            legend(ax(1), string(n_correct))
            plot(ax(2),t_subjective,Value_incorrect')
            legend(ax(2), string(n_incorrect))
            plot(ax(3),t_subjective,Value_bail')
            legend(ax(3),string(n_bail))
            plot(ax(4),t_subjective,RPE_correct')
            plot(ax(5),t_subjective,RPE_incorrect')
            plot(ax(6),t_subjective,RPE_bail')
            xlim(ax(4),get(ax(1),'xlim'))
            xlim(ax(5),get(ax(2),'xlim'))
            xlim(ax(6),get(ax(3),'xlim'))
        end
    end
end