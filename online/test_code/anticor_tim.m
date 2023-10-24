function [b r] = anticor_tim(varargin) 
% ANTICOR/TECHNICAL Modified Borodin-EL-Yaniv-Gogan ANTICOR algorithm 
%
% [B I] = ANTICOR(P,W,NORM,T,B0,DIR,PV,RS) Find the controls B for the next 
%   trading increment given the window size W, the index of the last 
%   trading increment T, historical sequence of market prices P and the 
%   current controls B0 (by the end of trading increment T). Matrix X has 
%   time as rows and stocks as columns i.e. Time x Entity. There are M 
%   stocks. NORM is by default TRUE to normalize weights to ONE, if FALSE 
%   weights are normalised to zero. DIR flag is set to +1 by default this 
%   can be set to -1, this reverses the direction of bet transfer between 
%   stock bets. DIR>0 implies the algorithm will exploit mean-reversal, 
%   DIR<0 implies the algorithm will exploit trend-following. PV Flag is 
%   TRUE to present-value the controls before the algorithm implementation, 
%   this is by default TRUE. RS is the returns to be used for attribution 
%   if this does not conform with the return of prices P. Typically this 
%   is empty or ignored, if P is smoothed then for attriubtion purposes RS 
%   will be the raw returns. I is the price index using the price P if no
%   RS is provided else RS is used.
%
% Note 1: NaN and missing data are ignored at the level of the algorithm
% Note 2: Long-only vs. Long-short is set using NORM
% Note 3: Direction of transfer can be changed using DIR
% Note 4: The control at time T, B(T), is the control at the beginning of
%   the T-th time period. The return R(T) is the return from the beginning of
%   the T-th time period until the end of the T-th time period. This ensures
%   that B(T) means that B(T) aligns with the prices at time T-1 but with the
%   return measured at time T. B(T) is invested at time T-1 and attributed
%   at time T as it is determined at time T-1 for the T-th period. I.E. 
%
%   R_B(T) = B(T)^T * R(T) =  B(T)^T * (P(T)/P(T-1))
%
% Example 1: Long-only at a single date-time
%     >> p = [ones(1,3); cumprod(exp(0.03*randn(10,3) + 0.01))];
%     >> [b1] = anticor(p,5,true)
%
% Example 2: Long-only sequence with initial BH conditions
%     >> p = [ones(1,3); cumprod(exp(0.03*randn(10,3) + 0.01))];
%     >> [b1] = anticor(p,5,true)
%
% Example 3: Long-short
%     >> p = [ones(1,3); cumprod(exp(0.03*randn(10,3) + 0.01))];
%     >> [b1] = anticor(p,w,false)
%     
%   See also: 

%   Author: Tim Gebbie, Raphael Nkomo, QT Capital Management

%   Reference: Bordin, A., El-Yaniv, R., Gogan, V., "Can we learn to beat 
%   the best stock", Journal of Artificial Intelligence Research 21 (2004) 
%   579-594

% TODO  : Modify to cope with missing data
% TODO  : Modify to cope with long-short
% FIXME : deal with repeated prices (row of zero returns)

% /CVSRepo/MATLAB/toolboxes/technical/functions/anticor.m,v 1.4 2008/07/11 10:56:42 Tim Gebbie Exp

%% Required data inputs
p = varargin{1}; % window size
w = varargin{2}; % last trading time-date
[n,m] = size(p);
tn=1:n-1;    % the time index
norm = true; % long-only (true) and long-short (false)
dir = +1;    % direction flag (+1 mean-reversals -1 for momentum)
type = 'Simple';
bias = []; % default bias is zero
leverage = 1; % default leverage is one

%% 0. Compute relative log-prices (returns)
x = tick2ret(p, [], type);
% repeated returns
ridx=transpose((sum(transpose((x==0)))==size(x,2)));
% compounding
switch type
    case 'Simple'
        x = (1+x); 
        %x=x;
    case 'Continuous'
        x = exp(x); % was exp(x)
    otherwise
end
% remove NaN
x(isnan(x)) = NaN;
% remove Inf
x(isinf(x)) = NaN;
% remove repeats
% FIXME x(ridx,:) = NaN;
rS = x; % attribution returns
%% Legacy control flags
pv_flag = true; % TRUE to present-value controls before ALGO

%% Input arguments
switch nargin
    case 2    
    case 3
        if islogical(varargin{3}), norm = varargin{3}; end;
        
    case 4
        if ~isempty(varargin{3}) & islogical(varargin{3}), norm = varargin{3}; end;
        tn = varargin{4}; % market price sequence
        [n,m] = size(p);
        if isempty(varargin{4}), tn=1:n-1; end;
        
    case 5
        if ~isempty(varargin{3}) & islogical(varargin{3}), norm = varargin{3}; end;
        % market price sequence
        if ~isempty(varargin{4}),  tn = varargin{4}; end;
        if length(tn)>1, error('T is expected to be of size 1 for online execution'); end;
        if ~isempty(varargin{5})
            b(tn,:) = varargin{5}; % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = nansum(b(tn,:)); % get the bias
            leverage = nansum(abs(b(tn,:))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial 
            % ---------------------------------------------------
        end
        
    case 6
        if ~isempty(varargin{3}) & islogical(varargin{3}), norm = varargin{3}; end;
        if ~isempty(varargin{4}),  tn = varargin{4}; end;  % market price sequence
        if ~isempty(varargin{5})
            b(tn,:) = varargin{5}; % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = nansum(b(tn,:)); % get the bias
            leverage = nansum(abs(b(tn,:))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial
            % ---------------------------------------------------
        end
        dir = varargin{6}; % direction
     
    case 7   
        if ~isempty(varargin{3}) & islogical(varargin{3}), norm = varargin{3}; end;
        if ~isempty(varargin{4}),  tn = varargin{4}; end;  % market price sequence
        if ~isempty(varargin{5}) & (size(varargin{5},1)==size(varargin{1},1)-1) % initial weight same size as price list
            b(tn,:) = varargin{5}(tn,:); % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = max(nansum(b(tn,:))); % get the bias
            leverage = max(nansum(abs(b(tn,:)))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial
            % ---------------------------------------------------
        elseif ~isempty(varargin{5})
            b(tn,:) = varargin{5}; % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = nansum(b(tn,:)); % get the bias
            leverage = nansum(abs(b(tn,:))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial
            % ---------------------------------------------------
        end
        if ~isempty(varargin{6}),  dir = varargin{6}; end; % direction
        if ~isempty(varargin{7}) & islogical(varargin{7}),  pv_flag = varargin{7}; end; % present-value flag
        
    case 8
        if ~isempty(varargin{3}) & islogical(varargin{3}), norm = varargin{3}; end;
        if ~isempty(varargin{4}),  tn = varargin{4}; end;  % market price sequence
        if ~isempty(varargin{5}) & (size(varargin{5},1)==size(varargin{1},1)-1) % initial weight same size as price list
            b(tn,:) = varargin{5}(tn,:); % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = max(nansum(b(tn,:))); % get the bias
            leverage = max(nansum(abs(b(tn,:)))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial
            % ---------------------------------------------------
        elseif ~isempty(varargin{5})
            b(tn,:) = varargin{5}; % current portfolio
            % ---------------------------------------------------
            % Initialise the bias and leverage
            bias = nansum(b(tn,:)); % get the bias
            leverage = nansum(abs(b(tn,:))); % get the leverage
            if leverage==0, leverage=1; end; % trap all zero weight initial
            % ---------------------------------------------------
        end
        if ~isempty(varargin{6}),  dir = varargin{6}; end; % direction
        if ~isempty(varargin{7}) & islogical(varargin{7}),  pv_flag = varargin{7}; end; % present-value flag
        if all(size(varargin{8})==size(x)), rS = varargin{8}; else error('rS dimensions incorrect'); end;
        
    otherwise
        error('Incorrect Input Arguments');
end

%% Long-only vs. Long-short initialisation
if nargin <5 | isempty(varargin{5}),
    % find missing data
    nidxn = ~isnan(x(1:max(2*w,tn(1)),:));
    % no preference initial conditions
    if norm
        % long-only has equal weights when there is no information
        for t=1:max(2*w,tn(1))
            % find the normalization
            nannorm = transpose(sum(nidxn(t,:)));
            % missing data initialize to zero
            b(t,~nidxn(t,:)) = 0;
            % non-missing data set to be fully invested
            b(t,nidxn(t,:)) = (1 ./ nannorm) * ones(1,nannorm); % equally weighted bets
        end
    else
        % long-short has not bets where there is no information
        b(1:max(2*w,tn(1)),:) = zeros(length(1:max(2*w,tn(1))),m); % no bets
    end
    % set the bias if it has not been set leverage is by default 1
    if isempty(bias)
        % set the bias if is has not already been set
        if norm
            bias = 1;
        else
            bias = 0;
        end
    end
end

%% Error check
% ensure the input weights are normalized
if abs(nansum(b(tn(1),:))-bias)>eps*100, error('Controls are incorrectly normalized'); end;

%% 1. Loop over time
for t=tn
    % 0. return the current portfolio if t < 2 w
    if t >= 2*w, % Inside the moving window

        % 1. Only present value if outside of the initial window
        if norm % Fully-Invested
            if pv_flag
                % present value fully-invested weights: weights are
                % renormalized because the portfolio is fully-invested
                % NaN weights are set to zero to ensure that NaN's only
                % arise from the return data directly
                b(t+1,:) = (b(t,:) .* x(t,:)) ./ nansum(b(t,:) .* x(t,:));
                b(t+1,isnan(b(t+1,:))) = 0; % ensure initial NaNs are ZERO
            else
                % no present value
                b(t+1,:) = b(t,:);
            end;
        elseif ~norm % Long-Short
            if pv_flag
                % present-value the weights: the weights are not reset to
                % zero, the bias nor leverage because it is long-short
                b(t+1,:) = b(t,:) .* x(t,:);
                b(t+1,isnan(b(t+1,:))) = 0; % ensure initial NaNs are ZERO
            else
                % no present value
                b(t+1,:) = b(t,:);
            end;
        end

        % 2.0 compute LX1(w,m), LX2(w,m), mu1 and mu2
        switch type
            case 'Simple'
                lx1 = x(t-2*w+1:t-w,:);
                lx2 = x(t-w+1:t,:); 
            case 'Continuous'
                lx1 = log(x(t-2*w+1:t-w,:));
                lx2 = log(x(t-w+1:t,:)); 
            otherwise
        end           
        
        % 2.1 compute the means and standard deviations
        mu1 = nanmean(lx1,1);
        mu2 = nanmean(lx2,1);
        if w>1
            s1  = nanstd(lx1,1);
            s2  = nanstd(lx2,1);
        else
            s1 = zeros(size(mu1));
            s2 = zeros(size(mu2));
        end
        
        % 2.2 remove the NaNs
        % remove single return time-series
        l1nidx = (sum(isnan(lx1))>0);
        l2nidx = (sum(isnan(lx2))>0);
        % remove NaN from means and standard deviations
        nidx1 = isnan(mu1);
        nidx2 = isnan(mu2);
        nidx3 = isnan(s1);
        nidx4 = isnan(s2);
        % logical index
        nidx  = ~(nidx1 | nidx2 | nidx3 | nidx4 | l1nidx | l2nidx);
        % numerical index to NaN's is reset
        nidxn = find(nidx);
        
        % 2.3 correct the data for NaN
        lx1 = lx1(:,nidx);
        lx2 = lx2(:,nidx);
        mu1 = mu1(nidx);
        mu2 = mu2(nidx);
        s1  = s1(nidx);
        s2  = s2(nidx);
        sum(nidx)

        % 3. compute the correlation matrix Mcor(i,j) and deal with the
        % case of window size w=1 for pure mean-reversals 
        if w>1
            Mcov = epsclean((1/(w-1)) * (lx1-(ones(w,1) * mu1))' * (lx2-(ones(w,1) * mu2)));
            Mcor =  Mcov ./ epsclean(s1' * s2); % Mcor(i,j) = Mcov(i,j) / s1(i) s2(j)
            Mcor(isinf(abs(Mcor)))=0; % remove zero sigma
            Mcor(isnan(abs(Mcor)))=0; % remove NaN data
            Mcor = epsclean(Mcor); % ensure small numbers are zero
        else
            Mcov = zeros(sum(nidx));
            Mcor = ones(sum(nidx));
        end
        
        % intialise claim(i to j) = 0
        claim = zeros(sum(nidx));
        % 4. compute claims for 1 <= i,j <= m
        for i=1:sum(nidx)
            for j=1:sum(nidx)
                % 5. if mu2(i) >= mu2(j) and Mcor(i,j)>0 then compute
                if (mu2(i)> mu2(j)) && (Mcor(i,j)>0) 
                    if w>1
                        claim(i,j) = Mcor(i,j) + max(-Mcor(i,i),0) + max(-Mcor(j,j),0); % claim i -> j
                    else
                        % pure mean reversion for w=1
                        claim(i,j) = Mcor(i,j) - Mcor(i,i) - Mcor(j,j); % claim i -> j
                    end
                end;
            end
        end
        % (a) claim(i to j) = claim(i to j) + Mcor(i,j)
        % (b) if Mcor(i,i)<0 then claim(i to j) = claim(i to j) - Mcor(i,i)
        % (c) if Mcor(j,j)<0 then claim(i to j) = claim(i to j) - Mcor(j,j)
        %
        % 6. compute the new portfolio
        % initialise controls (proportions) b(t+1,:) = b(t,:) for 1 <=i, j<=m
        transfer = zeros(sum(nidx));
        % compute the transfer function 
        for i=1:sum(nidx),
            if (sum(claim(i,:))~=0)
                for j=1:sum(nidx),
                    if (i~=j) && (sum(claim(i,:))~=0)
                        % (a) transfer(i to j) = b(t,i) * claim(i to j) / sum_j claim(i to j)
                        if norm
                            % sum of the transfers is 1
                            transfer(i,j) = dir * (b(t+1,nidxn(i)) * claim(i,j)) /  sum(claim(i,:));
                        else
                            % sum of transfers is zero maximum claim is 3
                            transfer(i,j) = 1/3 * dir * claim(i,j);
                        end
                    end
                end
            end
        end
        % 7. compute the new portfolio proportions
        for i=1:sum(nidx),
                % (b) b(t+1,i) = b(t+1,i) - transfer(i to j) + transfer(j to i) 
                b(t+1,nidxn(i)) = b(t+1,nidxn(i)) + (nansum(transfer(:,i)) - nansum(transfer(i,:))); % FIXME + <-> -
        end
        % renormalize weights to correct for missing data residuals
        if norm
            % 1. re-normalize to one
            b(t+1,:) = b(t+1,:) ./ nansum(b(t+1,:));
            % 3. re-bias if necessary
            b(t+1,:) = epsclean(bias * b(t+1,:));
        else
            % long and short index
            lnidx = (b(t+1,:)>0);
            snidx = (b(t+1,:)<0);
            % all none zero weights
            nznidx = lnidx | snidx;
            if sum(nznidx)==0
                b(t+1,:) = b(t,:);
            else
                % initial bias
                bias0 = nansum(b(t+1,:));
                % 1. re-normalize to zero
                b(t+1,nznidx) = b(t+1,nznidx) - (1/sum(nznidx)) * bias0;
                % 2. re-leverage
                b(t+1,:) = epsclean(leverage * b(t+1,:) ./ nansum(abs(b(t+1,:))));
                % 3. re-bias
                % b(t+1,lnidx) = b(t+1,lnidx) + (1/sum(lnidx)) * bias;
                % b(t+1,snidx) = b(t+1,snidx) - (1/sum(snidx)) * bias;
            end
        end
    end
end

%% post-process for output compatible with online functionality
% FIXME if ~isempty(find(ridx)), b(find(ridx)+1,:)=NaN; end;
if length(tn)==1, 
    b = b(tn+1,:); 
    r = NaN * b; 
else
    r = NaN * zeros(size(b));
    if ~isempty(rS), 
        % note that rS is a price relative => return = rS - 1
        xret = rS - 1;
        % must exclude the last control for which no return is avaliable
        % control b(n) at time n is invested at time n-1 for return
        % computed from n-1 until n. r(n) = sum_i b(n,i) * r(i,n)
        if size(b,2)==1
            % deal with trivial single stock case
            r = ret2tick(b(tn,:) .* xret(tn,:),[],[],[],type);
        else
            r = ret2tick(transpose(nansum(transpose(b(tn,:) .* xret(tn,:)))),[],[],[],type); 
        end
    end;
end



