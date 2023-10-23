function [b] = universal(x)
    % Universal algorithm for pairs of stocks
    % [b] = universal(x) finds portfolio weights b for the next day
    % given the sequence of price relatives x; row t of x consists of the
    % price relatives of the 2 stocks for day t
    %
    % Author: David White
    % T. M. Cover. Universal portfolios. Mathematical Finance, 1(1):1–29,
    % 1991.

    %dimensions of data tx2
    [t, m]=size(x);

    %1. initialise portfolio
    b=zeros(t,m);
    b(1,:)=1/m*ones(1, m);

    %2. for each day n
    for n=1:t-1
        num=0;
        denom=0;
        %2.a. for each i from 0,...,20
        for i=0:20
            %2.a.i. calculate Sn
            beta=[i/20 1-i/20]';
            Sn=prod(x(1:n,:)*beta);
            %2.a.ii. calculate numerator and denominator sums
            num=num+(i/20)*Sn;
            denom=denom+Sn;
        end
        %2.b. calculate portfolio for day n+1
        bn=num/denom;
        b(n+1,:)=[bn 1-bn];
    end

end