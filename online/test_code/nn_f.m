classdef nn_f
    % PATTERN Pattern matching and learning class
    %
    % The class implements both online and offline pattern matching and
    % learning over k-tuples for M objects and N features as specified in
    % an MxNxP data matrix X. The algorithm searches for L nearest-neighbour
    % K-tuple matches in the provided partition of data. A qaudratic
    % approximation is used to find the log-optimal portfolio using T+1
    % expected return subsequent to the pattern matching times T for each
    % k and ell to find H(K,L,CI;T) and SH(K,L,CI;T) for each K-tuple
    % L value and cluster CI. The controls are then aggregate using machine
    % learning to provide the controls B and realised accumulated
    % performance S. If the K matching patterns are predefined then K is
    % the index over the matching patterns. The default is and empty
    % matching pattern and to use the last k-tuple as the matching pattern.
    %
    % See Also: PATTERN/MATCH, PATTERN/OFFLINE, PATTERN/LEARN, PATTERN/ONLINE, HFTS
    
    %
    % A. OHLC patterns
    % B. Fundamental model patterns
    % C. (side information partitioning)
    %
    % 1. *Data* 
    %   2.1 M   (date-times)
    %   2.2 N  entities (objects e.g. stocks)
    %   2.3 P  features (OHLC)
    % 2. *Pattern* (k-tuple) [historic or user provided]
    %   2.1 k-tuple (k free parameter)
    %       2.2.1. k=1...K
    %       2.2.2. k=1,...,N*K multiples of N (DSI)
    %   2.2 k-tuple ConSet=[A,b] A*x>=b active/absolute (0,1)
    %       2.2.1. active (many stocks,stock+cash)
    %       2.2.2. absolute (single stock, long only portfolio)
    % 3. *Segment (s) and Partition (p)*
    %   3.1. [111...1] : full partition search for L nearest neighbours
    %   3.2. [10..0],[010...0],...,[0...01] : L partitions for single best fit in each
    %        partition
    %   3.3. [001],[011],[111] for partitions weight to current time
    %   3.4. ball radius (1-sigma ball relative to average distance)
    % 4. *Distance*
    %   4.1. surface norm
    %   4.2. vector norm
    % 5. *Predictors*
    %   5.1. dependent (correlated) matching times j(n,ell)
    %   5.2. independent matching times j(1...n,ell)
    %   5.3. different matching times ???
    %   5.4. [equi-probable] view is geometric average E[r]_t = exp(r_1+...+r_L)
    % 6. *Agents*
    %   6.1. Active/Absolute mean-variance
    %   6.2. Canonical Agents:
    %       6.2.1. controls H (NxMxT) M objects, N agents, T times)
    %       6.2.2. performance SH (arithmetic) (TxN for T times and N agents)
    %       6.2.3. horizon parameter (t to include delay)
    %   6.3. allow mixture of long-only with cash neutral.
    % 7. *Algorithm*
    %   7.1. Online/Offline
    %   7.2. Parallel computing
    %   7.3. Interface agents (h,SH)
    % 8. *Machine Learning* (EG,EW,EWMA,UNIV + abs/act + part/comb)
    %   8.1. Weighted Arithmetic Average over all agents to create optimal predictors
    %       8.1.1. Performance weighted averaging (using arithmetic averaging)
    %           B(T) = SUM(K,L) (SH(T-1|K,L) H(T|K,L)) / SUM(K,L) SH(T-1|K,L)
    %           this is probability wieghted where the probability is
    %           proportional to the returns.
    %       8.1.2. Exponentially weight more recent performance data
    %   8.2. parameters: window W_L, forgetting factor Lamba_L
    %   8.3. either fully invested or active.
    % 9. *Sectors and States* (clusters)
    %   9.1. CI PxN for P clusters and N objects
    
    % Author: Tim Gebbie
    
    %% public properties
    properties
        b = []; % aggregated controls    (T,M) Time x Objects
        S = []; % aggregated performance (T,1) Time x 1
        h = []; % agent controls      (N,M,T)  Agents x Objects x Time
        SH = []; % agents performance (T,N)    Time x Agents
        %{
    end
    %% private properties
    properties (Access = private)
        %}
        x = []; % data as price relatives (M,F,T) Objects x Features x Time
        k = []; % k-tuple size
        ell = []; % partition size
        mnp = []; % dimensionality (M,F,T)
        qH = []; % agents controls (T,N) Time x Agents
    end
    %% methods
    methods
        %% constructor
        function p = nn_f(x, k, ell)
            % P = PATTERN/PATTERN constructor
            %
            % P = PATTERN(X,K,ELL,CI,NTYPE,LNTYPE,LTYPE) For data X (price relatives)
            %   and is MxNxP dimensional data for M objects, N features and
            %   P date-times. dim(X) is dim(Price)-1. A typical object is
            %   a stocks, e.g. AGL, a typical feature is a factors such as
            %   OPEN, HIGH, LOW, CLOSE (OHLC), and date-times are the time
            %   stamps for the data.
            %
            % Example 1:
            % >> x = 1+0.2*randn(10,1,1000);
            % >> p = offline(p);
            % >> p.S
            %
            %  P = PATTERN(X,XNK,ELL,CI,NTYPE,LNTYPE,LTYPE) For matching
            %   pattern XNK instead of K-tuple range.
            %
            %  P = PATTERN(X,XNK,ELL,CI,NTYPE,LNTYPE,LTYPE,TREND)
            %   LONGTERM is TRUE to use long-term relative views. This is 
            %   by default false.
            %
            % The data X is homogenize in time. K is the set of tuple
            % sizes if there is no matching pattern e.g. [1:3], if there
            % is a matching pattern it is the number of matching patterns. 
            % The number of matching times is ELL this is typically 10 
            % per partition. The cluster definitions is CI, this is by 
            % default the trivial cluster (all the objects) when CI is 
            % empty. CI is KxM for K clusters of M objects. NTYPE is the 
            % strategy normalisation and can be either 'active', 'absolute' 
            % or 'gyorfi_opt'. LNTYPE is the agent normalisations and can
            % be either 'active' or 'absolute'.
            % The LTYPE is the learning type this is by default 'univ'.
            %
            % Note: X are price relatives (P(t)/P(t-1)). These can be
            % conveniently computed using EXP(DIFF(LOG(P))).
            %
            % See Also PATTERN/DISPLAY, PATTERN/SUBSREF
            
                %% Compute price relatives as returns
                p.x=x;
                p.k=k;
                p.ell=ell;
                p.mnp = size(p.x);
                % control parameters
                L = size(p.ell,2); % partitions (param 1)
                K = size(p.k,2); % tuples (param 2)
                % initialise state variables and pre-allocate memory
                p.b = nan(p.mnp(3)+1,p.mnp(1)); % Time x Objects        
                p.h = nan(L*K,p.mnp(1),p.mnp(3)+1); % Agents x Objects x Time    
                p.qH = nan(p.mnp(3)+1,L*K); % Time x Agents
                p.S = ones(p.mnp(3)+1,1); % Time x 1
                p.SH = ones(p.mnp(3)+1,L*K); % Time x  Agents
        end
        %% offline (to force offline estimation)
        function p=offline(p)
            % PATTERN/OFFLINE Offline estimation
            %
            % P = OFFLINE(P) to estimate (H,SH) for T=K*L:T0 for K and L.
            %
            % See Also: PATTERN/ONLINE
            
            tmin0 = 3*max(p.ell);
            t0 = p.mnp(3);
            if t0<tmin0
                error('pattern:offline','Not enough Data L*K>T');
            end
            % Find matching times j for pattern [offline loop]
            for t=tmin0:t0 % time loop
                %% control parameters
                L = size(p.ell,2); % partitions (param 1)
                K = size(p.k,2); % tuples (param 2)
                %% Model Identification Loop
                % initial controls
                % t+1 rows as the last row is for an unrealised return
                ht = nan(K*L,p.mnp(1)); % agents controls
                xt = p.x(:,:,1:t);
                % exit if there is not enough data
                if 3*max(p.ell)>t
                    error('pattern:online','Not enough data K*L>T');
                end
                if (t>p.mnp(3))
                    error('pattern:online','Not enough data T> dim(X,T)');
                end
                % matching loop
                xtw0 = xt;
                for ell0=1:L % parameter 1 -> passed to match (ell neighbours)
                    for k0=1:K  % parameter 2 -> passed to match (k-tuple)
                        % expert index KLrow(w0,k0,ell0)
                        KLrow = sub2ind([K,L],K-k0+1,ell0);
                        % (k,ell)-agents matched pattern for cluster ci(w0)
                        % --- 1. select the pattern -------------------------
                        % k index tuple size
                        xnk0 = xtw0(:,:,end-p.k(k0)+1:end)-1; % r= R - 1
                        % --- 2. pattern matching ---------------------------
                        %% Loop parameters wrt to partition
                        m0 = size(xtw0,1); % number of objects
                        % the pattern as returns computed from price relatives
                        k0 = size(xnk0,3);
                        % check consistency
                        if any(size(k0)>1), error('portchoice:pattern:match','Incorrect K'); end
                        if any(size(ell0)>1), error('portchoice:pattern:match','Incorrect L'); end

                        %% Find matching times j for pattern
                        % the partition
                        psi = true(1,t);
                        % for agent h(k,ell) for cluster ci, and partition psi
                        snk = xtw0(:,:,psi)-1;
                        % find the times (this allows for inhomogenous partitions)
                        jell = find(psi);
                        % reset the distance measure for partition
                        ed = Inf(size(snk,3),m0);
                        %% get the test tuples by looping over partition snk
                        for j=k0:size(snk,3) - 1 % only allow jell+1 in partition
                            % distance (element-wise difference)
                            edi = xnk0 - snk(:,:,j-k0+1:j);

                            if (k0==1)
                                % 2-norm computed sum_i sqrt(a_i^2)
                                ed(j,1:m0) = norm(edi);
                            else
                                % reshape the distance by objects and factors (Edited by FL)
                                % computes distance like Gyorfi
                                edi = reshape(edi,size(edi,1),size(edi,2)*size(edi,3));
                                ed(j,1:m0) = norm(edi);
                                % 2-norm computed for each object independently (Edited by FL)
                                %             for a=1:m0
                                %                 ed(j,a) = norm(edi(a,:))
                                %                 norm(edi(a,:))
                                %                 norm(edi)
                                %             end
                            end
                        end % j segment loop
                        %% sort the matching times

                        % Edited by FL

                        pl = 0.02 + 0.5*((ell0-1)/(L-1));
                        pell0 = floor(pl*j);


                        % ell matches in a single partition
                        [~,nj]=sort(ed,'ascend'); % ~ -> pat.ed not required

                        % Take into account ties of the norm (Edited by FL)
                        if pell0>0
                            while pell0<length(nj) && (ed(nj(pell0))==ed(nj(pell0+1)))
                                pell0 = pell0 + 1;
                            end
                        end

                        % update the matching times (Edited by FL)
                        njell = jell(nj(1:pell0,:));
                        % update the norm (Edited by FL)
                        ed = ed(nj(1:pell0,:));

                        % update the distances
                        % pat.ed = pat.ed(nj(pat.ell0,:));
                        % find the ell matching times
                        % --- LOOK AHEAD RULE ------------
                        jnk = njell+1; % 1-step look ahead
                        %% find the predictions (Edited by FL)
                        E = eye(m0,m0);
                        hatx = zeros(pell0,m0);
                        % initialise prediction vector
                        for mj = 1:m0 % loop over objects
                            % select first factor (price relative) @JNK matching times
                            hatx(:,mj) = xtw0(mj,1,jnk(:,mj))-1; % r = R - 1
                            % matching error (diag covariance matrix)
                            E(mj,mj) = mean(ed(:,mj).^2);
                        end % loop over objects/stocks
                        xx = hatx +1;
                        [~, m] = size(xx);
                        b0 = (1/m)*ones(1,m);
                        optfun = @(b)(-prod(xx*transpose(b)));
                        [hkl2, ~] = fmincon(optfun,b0,[-eye(m);eye(m)],...
                            [zeros(1,m)';ones(1,m)'],ones(1,m),1,[],[],[],...
                            optimset('Algorithm','sqp','Display','off'));
                        hklt = hkl2';
                        % ------------------------------------------------
                        % expert controls per cluster mapping
                        ht(KLrow,:) = transpose(hklt);
                    end % k
                end % ell
                % initialise h,SH if it is the trivial object
                % compute the update performance for the prior agent step
                if any(isnan(p.h(:,:,t)))
                    dSH = ones(size(p.SH(t,:)));
                else
                    dSH = transpose((p.h(:,:,t) * (p.x(:,1,t)-1))+1); % was exp
                end
                % remove NaN
                ht(isnan(ht))=0;
                % update the agents
                p.h(:,:,t+1) = ht;
                % update the agent accumulate performance (geometric returns)
                p.SH(t+1,:) =  p.SH(t,:) .* dSH;
                %% online update the learning
                %% Machine Learning
                if (size(p.h,1)==1)
                    % only single agent
                    % ---- ONLINE update ----
                    b0 = p.h(1,:,t+1);
                    qH0 = 1;
                    % -----------------------
                else
                    % multiple agents
                    qH0 = p.SH(t+1,:);
                    %% renormalise the weights
                    qH0 = qH0 ./ sum(qH0);
                    %% create the performance weighted combination of experts.
                    % ONLINE
                    % -------------------------------------------------------
                    b0 = qH0 * p.h(:,:,t+1);
                    % -------------------------------------------------------
                    %% compute normalization abs(long) + abs(short)
                    tb = nansum(abs(b0));
                    % renormalize controls (leverage=1) [FIXME LEV]
                    if tb==1
                    elseif tb>eps
                        b0 = (1/tb) * b0;
                        % update the agent mixture weights for leverage
                        qH0 = (1/tb) * qH0;
                    else
                        % update the agent mixture weights for leverage
                        qH0 = zeros(size(qH0)); % FIXME (should be equally weighted)
                        % zero weights
                        b0 = zeros(size(b0));
                    end
                end % only 1 agent
                %% compute the leverage corrected output price relative
                % compute the updated returns
                dSH = p.SH(t+1,:)./p.SH(t,:);
                % uses the inputed online structure reference.
                if all(isnan(p.qH(t,:)))
                    dS = 1;
                else
                    dS = ((dSH-1) * transpose(p.qH(t,:)))+1; %  LINRET was exp
                end
                % update the properties
                p.qH(t+1,:) = qH0;
                p.b(t+1,:) = b0;
                p.S(t+1,:) = p.S(t,:) * dS;
            end
            
        end
        %% online (to force online estimation)
        %{
        function p=online(p, t)
            % PATTERN/ONLINE Offline estimation
            %
            % P = ONLINE(P) to estimate (H,SH,B,S) at T for the range of
            %   K and L over the specified clusters CI. This requires the
            %   online structure for learning to have been initialised. T
            %   is taken to be the last time in the object.
            %
            % P = ONLINE(P,T) to estimate online values at time T using the
            %   data from times 1 to T (1:T).
            %
            % See Also: PATTERN/ONLINE, PATTERN/MATCH, PATTERN/MATCH
            %% control parameters
            L = size(p.ell,2); % partitions (param 1)
            K = size(p.k,2); % tuples (param 2)
            %% Model Identification Loop
            % initial controls
            % t+1 rows as the last row is for an unrealised return
            ht = nan(K*L,p.mnp(1)); % agents controls
            xt = p.x(:,:,1:t);
            % exit if there is not enough data
            if 3*max(p.ell)>t
                error('pattern:online','Not enough data K*L>T');
            end
            if (t>p.mnp(3))
                error('pattern:online','Not enough data T> dim(X,T)');
            end
            % matching loop
            xtw0 = xt;
            for ell0=1:L % parameter 1 -> passed to match (ell neighbours)
                for k0=1:K  % parameter 2 -> passed to match (k-tuple)
                    % expert index KLrow(w0,k0,ell0)
                    KLrow = sub2ind([K,L],K-k0+1,ell0);
                    % (k,ell)-agents matched pattern for cluster ci(w0)
                    % --- 1. select the pattern -------------------------
                    % k index tuple size
                    xnk0 = xtw0(:,:,end-p.k(k0)+1:end)-1; % r= R - 1
                    % --- 2. pattern matching ---------------------------
                    [hklt]=match(xtw0,t,xnk0,ell0,L); % online
                    % ------------------------------------------------
                    % expert controls per cluster mapping
                    ht(KLrow,:) = transpose(hklt);
                end % k
            end % ell
            % initialise h,SH if it is the trivial object
            % compute the update performance for the prior agent step
            if any(isnan(p.h(:,:,t)))
                dSH = ones(size(p.SH(t,:)));
            else
                dSH = transpose((p.h(:,:,t) * (p.x(:,1,t)-1))+1); % was exp
            end
            % remove NaN
            ht(isnan(ht))=0;
            % update the agents
            p.h(:,:,t+1) = ht;
            % update the agent accumulate performance (geometric returns)
            p.SH(t+1,:) =  p.SH(t,:) .* dSH;
            %% online update the learning
            %% Machine Learning  
            if (size(p.h,1)==1)
                % only single agent
                % ---- ONLINE update ----
                b0 = p.h(1,:,t+1);
                qH0 = 1;
                % -----------------------    
            else
                % multiple agents
                qH0 = p.SH(t+1,:);
                %% renormalise the weights
                qH0 = qH0 ./ sum(qH0);
                %% create the performance weighted combination of experts.
                % ONLINE
                % -------------------------------------------------------
                b0 = qH0 * p.h(:,:,t+1);
                % -------------------------------------------------------
                %% compute normalization abs(long) + abs(short)
                tb = nansum(abs(b0));
                % renormalize controls (leverage=1) [FIXME LEV]
                if tb==1
                elseif tb>eps
                    b0 = (1/tb) * b0;
                    % update the agent mixture weights for leverage
                    qH0 = (1/tb) * qH0;
                else
                    % update the agent mixture weights for leverage
                    qH0 = zeros(size(qH0)); % FIXME (should be equally weighted)
                    % zero weights
                    b0 = zeros(size(b0));
                end
            end % only 1 agent
            %% compute the leverage corrected output price relative
            % compute the updated returns
            dSH = p.SH(t+1,:)./p.SH(t,:);
            % uses the inputed online structure reference.
            if all(isnan(p.qH(t,:)))
                dS = 1;
            else
                dS = ((dSH-1) * transpose(p.qH(t,:)))+1; %  LINRET was exp
            end
            % update the properties
            p.qH(t+1,:) = qH0;
            p.b(t+1,:) = b0;
            p.S(t+1,:) = p.S(t,:) * dS;
        end
        %}
        %% learning
        %{
        function p=learn(p, t)
            % PATTERN/LEARN Machine Learning based on performance
            %
            % P = LEARN(P) The updates the aggregated agents and agent
            %   performance (B,S) using the specified learning type.
            %
            % P = LEARN(P,TYPE) will reset the learning
            %
            % Table 1: Learning types
            % +----------------+---------------------------------------+
            % | TYPE           | Description                           |
            % +----------------+---------------------------------------+
            % | 'univlearning' | Universal learning (log-optimal)      |
            % |                | [PARAM=[NONE]                         |
            % | 'eglearning'   | Exponentiated Gradient                |
            % |                | [PARAM=ETA in [0,1] typ. [0,0.2]      |
            % | 'ewlearning'   | EWMA in control based learning        |
            % |                | [PARAM=LAMBDA in [0,1] typ. [0.9,0.99]|
            % +----------------+---------------------------------------+
            %
            % References:
            % [1] Cover, T., M. (1991) Universal Portfolios
            % [2] Gyorfi, L., Udina, F., Walk, H., (2008) Experiments on universal
            %           portfolio selection using data from real markets
            % [3] Cover, T. M., (1996), Universal Portfolios with Side Information
            % [4] Algoet, P. H., Cover, T. M., (1980) Asymptotic optimality and
            %           symptotic equipartition properties of log-optimum investments
            % [5] Helmbold, D., P., Schapire, R., E., Singer, Y., Warmuth, M.,
            %           K.,(1998) On-line portfolio selection using multiplicative updates
            %
            % See Also: PATTERN/MATCH
            
            % SHX = exp(diff(log(SH)) for price path SH!
            %% Machine Learning  
            if (size(p.h,1)==1)
                % only single agent
                % ---- ONLINE update ----
                b0 = p.h(1,:,t+1);
                qH0 = 1;
                % -----------------------    
            else
                % multiple agents
                qH0 = p.SH(t+1,:);
                %% renormalise the weights
                qH0 = qH0 ./ sum(qH0);
                %% create the performance weighted combination of experts.
                % ONLINE
                % -------------------------------------------------------
                b0 = qH0 * p.h(:,:,t+1);
                % -------------------------------------------------------
                %% compute normalization abs(long) + abs(short)
                tb = nansum(abs(b0));
                % renormalize controls (leverage=1) [FIXME LEV]
                if tb==1
                elseif tb>eps
                    b0 = (1/tb) * b0;
                    % update the agent mixture weights for leverage
                    qH0 = (1/tb) * qH0;
                else
                    % update the agent mixture weights for leverage
                    qH0 = zeros(size(qH0)); % FIXME (should be equally weighted)
                    % zero weights
                    b0 = zeros(size(b0));
                end
            end % only 1 agent
            %% compute the leverage corrected output price relative
            % compute the updated returns
            dSH = p.SH(t+1,:)./p.SH(t,:);
            % uses the inputed online structure reference.
            if all(isnan(p.qH(t,:)))
                dS = 1;
            else
                dS = ((dSH-1) * transpose(p.qH(t,:)))+1; %  LINRET was exp
            end
            % update the properties
            p.qH(t+1,:) = qH0;
            p.b(t+1,:) = b0;
            p.S(t+1,:) = p.S(t,:) * dS;
        end
        %}
    end % end methods
end
%% helper functions
%{
function [hklt] = match(xtw0, t, xnk0, ell0, L)
% PATTERN/MATCH Pattern Match for a given K-tuple and matching data set.
%
% [HKL]=MATCH(X0,P0,XNK,L0,NTYPE,TREND,HORZ) X0 are the factor relatives. P0 the 
%   partition. K0 is the k-tuple size over the current data set, it is 
%   computed from the user definied pattern (k-tupel) XNK. XNK is a MxNxK0 
%   double for X0 a MxNxP size double for T>>K0. L0 is the number of 
%   neighbours to include from the partitioning P0 of the data. NTYPE is 
%   the normalisation type. This excludes the time loop over the the 
%   matching horizon. This is the online version of the pattern matching 
%   and learning algorithm. HKL is a Mx1 control vector that satifies the 
%   normalisation type NTYPE. 
%
% See Also: PATTERN/ONLINE, PATTERN/OFFLINE, QUADBET

% Author: Tim Gebbie

%% Loop parameters wrt to partition
m0 = size(xtw0,1); % number of objects
% the pattern as returns computed from price relatives
k0 = size(xnk0,3);
% check consistency
if any(size(k0)>1), error('portchoice:pattern:match','Incorrect K'); end
if any(size(ell0)>1), error('portchoice:pattern:match','Incorrect L'); end

%% Find matching times j for pattern
    % the partition
    psi = true(1,t);
    % for agent h(k,ell) for cluster ci, and partition psi
    snk = xtw0(:,:,psi)-1; 
    % find the times (this allows for inhomogenous partitions)
    jell = find(psi);
    % reset the distance measure for partition
    ed = Inf(size(snk,3),m0);
    %% get the test tuples by looping over partition snk
    for j=k0:size(snk,3) - 1 % only allow jell+1 in partition
        % distance (element-wise difference)
        edi = xnk0 - snk(:,:,j-k0+1:j);
        
        if (k0==1)
            % 2-norm computed sum_i sqrt(a_i^2)
            ed(j,1:m0) = norm(edi);
        else
            % reshape the distance by objects and factors (Edited by FL)
            % computes distance like Gyorfi
            edi = reshape(edi,size(edi,1),size(edi,2)*size(edi,3));
            ed(j,1:m0) = norm(edi);
            % 2-norm computed for each object independently (Edited by FL)
%             for a=1:m0
%                 ed(j,a) = norm(edi(a,:))
%                 norm(edi(a,:))
%                 norm(edi)
%             end
        end
    end % j segment loop
    %% sort the matching times
    
    % Edited by FL
    
    pl = 0.02 + 0.5*((ell0-1)/(L-1));
    pell0 = floor(pl*j);
   
    
    % ell matches in a single partition
    [~,nj]=sort(ed,'ascend'); % ~ -> pat.ed not required

    % Take into account ties of the norm (Edited by FL)
    if pell0>0
        while pell0<length(nj) && (ed(nj(pell0))==ed(nj(pell0+1)))
            pell0 = pell0 + 1;
        end
    end

    % update the matching times (Edited by FL)
    njell = jell(nj(1:pell0,:));
    % update the norm (Edited by FL)
    ed = ed(nj(1:pell0,:));

    % update the distances
    % pat.ed = pat.ed(nj(pat.ell0,:));
    % find the ell matching times
    % --- LOOK AHEAD RULE ------------
    jnk = njell+1; % 1-step look ahead
%% find the predictions (Edited by FL)
E = eye(m0,m0);
hatx = zeros(pell0,m0);
% initialise prediction vector
for mj = 1:m0 % loop over objects
    % select first factor (price relative) @JNK matching times
    hatx(:,mj) = xtw0(mj,1,jnk(:,mj))-1; % r = R - 1
    % matching error (diag covariance matrix)
    E(mj,mj) = mean(ed(:,mj).^2);
end % loop over objects/stocks
xx = hatx +1;
[~, m] = size(xx);
b0 = (1/m)*ones(1,m);
optfun = @(b)(-prod(xx*transpose(b)));
[hkl2, ~] = fmincon(optfun,b0,[-eye(m);eye(m)],...
    [zeros(1,m)';ones(1,m)'],ones(1,m),1,[],[],[],...
    optimset('Algorithm','sqp','Display','off'));
hklt = hkl2';
end
%}