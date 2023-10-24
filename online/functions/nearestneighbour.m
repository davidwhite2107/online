function [b, S, h, SH] = nearestneighbour(x, k, ell)
    % Gyorfi-Udina-Walk's nearest neighbour algorithm
    % [b, S, h, SH]=nearestneighbour(x, k, ell) finds portfolio weights b for the next day given
    % the sequence of k values and the sequence of ell values that define
    % the agents and the sequence of price relatives x; row t of 
    % x consists of the price relatives of the m stocks for day t. S is the
    % performance of the portfolio over all the days, h consists of the
    % portfolio weights of each agent for each day and SH is each agents
    % performance for each day.
    %
    %
    % Author: David White
    % Reference: F. Udina & H. Walk L. Gyorfi. Nonparametric nearest 
    % neighbor based empirical portfolio selection strategies. 
    % Statistics & Risk Modeling, 26(2):145–157, 2008

    %% 1. initialise variables
    %dimensions of x
    [t0, m0]=size(x);

    %hyperparameters
    %number of matches
    L = size(ell,2);
    %length of window
    K = size(k,2);

    % initialise variables to be returned
    b = nan(t0+1,m0); % time x stocks
    h = nan(L*K,m0,t0+1); % agents x stocks x time
    qH = nan(t0+1,L*K); % time x agents
    S = ones(t0+1,1); % time x 1
    SH = ones(t0+1,L*K); % time x  agents

    % ensure there is enough data
    tmin0 = 3*max(ell);
    if t0<tmin0
        error('Error: Not enough Data.');
    end
    % first 3L portfolios are uniform
    b(1:tmin0,:)=(1/m0)*ones(tmin0,m0);
    % performance for first 3L days
    port_ret=zeros(tmin0,1);
    for t=1:tmin0
        port_ret(t)=b(t,:)*x(t,:)';
    end
    S(1:tmin0,:)=cumprod(port_ret);
    %% 2. loop over t
    for t=tmin0:t0 % time loop
        % initialise agents controls for day t
        ht = nan(K*L,m0);
        % data available at time t
        xt = x(1:t,:);

        for ell0=1:L % loop over ell
            for k0=1:K  % loop over k
                %% 3. find ell_hat nearest neighbours
                % 3.a. calculate euclidean norms of all possible matches
                % index for each agent's portfolio in ht
                KLrow = sub2ind([K,L],K-k0+1,ell0);
                % pattern to be matched
                xnk0 = xt(end-k(k0)+1:end,:)-1;
                snk = xt-1;
                jell = 1:t;
                % initialise distances
                ed = Inf(size(xt,1),m0);

                % calculate distances (Euclidean norm)
                for j=k0:size(xt,1) - 1 % only allow jell+1 in partition
                    % distance (element-wise difference)
                    edi = xnk0 - snk(j-k0+1:j,:);

                    if (k0==1)
                        % euclidean norm
                        ed(j,1:m0) = norm(edi);
                    else
                        % reshape edi so euclidean norm can be used
                        edi = reshape(edi,size(edi,1),size(edi,2));
                        ed(j,1:m0) = norm(edi);
                    end
                end

                % 3.b. get only the ell_hat nearest neighbours
                % calculate ell_hat
                pl = 0.02 + 0.5*((ell0-1)/(L-1));
                pell0 = floor(pl*j);
                % ell_hat nearest neighbour indices
                [~,nj]=sort(ed,'ascend');
                % take ties into account
                if pell0>0
                    while pell0<length(nj) && (ed(nj(pell0))==ed(nj(pell0+1)))
                        pell0 = pell0 + 1;
                    end
                end
                % matching times
                njell = jell(nj(1:pell0,:));
                jnk = njell+1;
                %% 4. find optimal portfolio based on nearest neighbours
                hatx = zeros(pell0,m0);
                for mj = 1:m0 
                    hatx(:,mj) = xt(jnk(:,mj),mj);
                end 
                %initial guess
                b0 = (1/m0)*ones(1,m0);
                %optimise
                optfun = @(b)(-prod(hatx*transpose(b)));
                [hkl2, ~] = fmincon(optfun,b0,[-eye(m0);eye(m0)],...
                    [zeros(1,m0)';ones(1,m0)'],ones(1,m0),1,[],[],[],...
                    optimset('Algorithm','sqp','Display','off'));
                ht(KLrow,:) = hkl2;
            end % k
        end % ell
        
        %% 5. Combine agents into one portfolio
        if any(isnan(h(:,:,t)))
            dSH = ones(size(SH(t,:)));
        else
            dSH = transpose((h(:,:,t) * (x(t,:)'-1))+1);
        end
        % remove NaN
        ht(isnan(ht))=0;
        % update the agents 
        h(:,:,t+1) = ht;
        % update the agent cumulative returns
        SH(t+1,:) =  SH(t,:) .* dSH;

        if (size(h,1)==1)
            % only single agent
            b0 = h(1,:,t+1);
            qH0 = 1;
        else
            % multiple agents
            qH0 = SH(t+1,:);
            % renormalise the weights
            qH0 = qH0 ./ sum(qH0);
            % create the performance weighted combination of experts.
            b0 = qH0 * h(:,:,t+1);
            % compute normalization
            tb = nansum(abs(b0));
            % renormalize weights
            if tb==1
            elseif tb>eps
                b0 = (1/tb) * b0;
                qH0 = (1/tb) * qH0;
            else
                % update the agent mixture weights for leverage
                qH0 = zeros(size(qH0));
                % zero weights
                b0 = zeros(size(b0));
            end
        end
        % compute the updated returns
        dSH = SH(t+1,:)./SH(t,:);
        if all(isnan(qH(t,:)))
            dS = 1;
        else
            dS = ((dSH-1) * transpose(qH(t,:)))+1; 
        end
        % update the properties
        qH(t+1,:) = qH0;
        b(t+1,:) = b0;
        S(t+1,:) = S(t,:) * dS;
    end
end