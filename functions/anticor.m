function [b]=anticor(x, w)
    % ANTICOR algorithm
    % [b]=anticor(x, w) finds portfolio weights b for the next day given
    % the window size w and the sequence of price relatives x; row t of 
    % x consists of the price relatives of the m stocks for day t
    %
    % Author: David White
    % Reference: Bordin, A., El-Yaniv, R., Gogan, V., "Can we learn to beat
    % the best stock", Journal of Artificial Intelligence Research 21 (2004)
    % 579-594


    %dimensions of data
    [tn,m]=size(x);

    %initialise portfolio
    b=zeros(tn,m);
    b(1, :)=1/m*ones(m,1);
    for t=1:(tn-1)
        
        %1. if t<2w use 1/m portfolio
        if t<2*w
            b(t+1, :)=b(t,:);
        else
            %2. initialise b_(t+1) to bhat_t
            b(t+1,:) = (b(t,:) .* x(t,:)) ./ nansum(b(t,:) .* x(t,:));

            %3.a. calculate LX1, LX2, mu1 and mu2
            lx1 = log(x(t-2*w+1:t-w,:));
            lx2 = log(x(t-w+1:t,:));
            mu1 = mean(lx1,1);
            mu2 = mean(lx2,1);

            %3.b. calculate sd1 and sd2
            if w>1
                s1  = std(lx1,1);
                s2  = std(lx2,1);
            else
                s1 = zeros(size(mu1));
                s2 = zeros(size(mu2));
            end
       
            %3.c. calculate Mcov and Mcor
            if w>1
                Mcov = epsclean((1/(w-1)) * (lx1-(ones(w,1) * mu1))' * (lx2-(ones(w,1) * mu2)));
                Mcor =  Mcov ./ epsclean(s1' * s2); % Mcor(i,j) = Mcov(i,j) / s1(i) s2(j)
                Mcor(isinf(abs(Mcor)))=0; % remove zero sigma
                Mcor(isnan(abs(Mcor)))=0; % remove NaN data
                Mcor = epsclean(Mcor); % ensure small numbers are zero
            else
                Mcov = zeros(m);
                Mcor = ones(m);
            end

            % 4.a. intialise claim_ij = 0
            claim = zeros(m);

            % 4.b. for i,j from 1:m, calculate claim_ij
            for i=1:m
                for j=1:m
                    if (mu2(i)> mu2(j)) && (Mcor(i,j)>0) 
                        if w>1
                            claim(i,j) = claim(i,j)+Mcor(i,j) + max(-Mcor(i,i),0) + max(-Mcor(j,j),0); % claim i -> j
                        else
                            % pure mean reversion for w=1
                            claim(i,j) = claim(i,j)+Mcor(i,j) - Mcor(i,i) - Mcor(j,j); % claim i -> j
                        end
                    end
                end
            end

            %5.a. initialise transfer_ij = 0
            transfer = zeros(m);

            %5.b. for i,j from 1:m, i != j, calculate transfer_ij
            for i=1:m
                if (sum(claim(i,:))~=0)
                    for j=1:m
                        if (i~=j) && (sum(claim(i,:))~=0)
                                transfer(i,j) = (b(t+1,i) * claim(i,j)) /  sum(claim(i,:));
                        end
                    end
                end
            end

            %6. compute portfolio weights for day t+1
            for i=1:m
                b(t+1,i) = b(t+1,i) + (sum(transfer(:,i)) - sum(transfer(i,:)));
            end

            %7. ensure portfolio is fully invested
            b(t+1,:) = b(t+1,:) ./ nansum(b(t+1,:));
            %convert small numbers to zero
            b(t+1,:)=epsclean(b(t+1,:));

        end
    end
end

