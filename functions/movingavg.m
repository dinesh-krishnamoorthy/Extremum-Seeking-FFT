function y = movingavg(x,n,order)
warning off
if nargin<3
    order = 0;
end
ne = numel(x);
n2 = round(n/2);
y = zeros(size(x));
for i=1:ne
    if i<n2 % before half of the window
        X  = x(1:n);
        ie = i; % Evaluation point
    elseif i>=(ne-n2) % less than half of the window left
        X  = x(end-n+1:end);
        ie = i-(ne-n); 
        %ie = i-(ne-n2); 
    else
        X  = x(i-n2+1:i+n2);
        ie = n2;
    end
    if order == 0 
        y(i) = mean(X);
    else
        p = polyfit([1:n]',X(:),order);
        y(i) = polyval(p,ie);
    end
end
warning on
    
    