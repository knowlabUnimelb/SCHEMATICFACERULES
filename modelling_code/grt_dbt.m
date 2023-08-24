% Compute drift rate (p(Step to A)) given dimension boundaries, means and
% variances
function [p] = grt_dbt(x, dx, sx, coactive)

% x  = [15    9.5  ]; % [dimx dimy]
% dx = [15.97 10.03]; % [bndx bndy]
% sx = [1.89  0.988]; % [sigmax sigmay]

if ~coactive
    if size(x, 1) > 1
        dx = repmat(dx, size(x, 1), 1);
        sx = repmat(sx, size(x, 1), 1);
    end
    
    zx = ((dx - x)./sx);
    p  = normcdf(zx, 0, 1);
    
    % Truncate drift rate
    p(p<.00001) = .00001;
    p(p>.99999) = .99999;
    % p2 = ztop(zx);

else
    %   B | D
    %  -------
    %   A | C

   % Check positive semi-definiteness 
   [~,err] = cholcov(diag(sx),0);
   if err~=0
       warning('Covariance matrix not positive definate. Likely cause is that the variances have grown unusually large.');
       save errortemp
   end
    
    aRegion = [-inf -inf; dx(1) dx(2)];
    bRegion = [-inf dx(2); dx(1)  inf];
    cRegion = [dx(1) -inf;  inf dx(2)];
    dRegion = [dx(1) dx(2);  inf  inf];
    
    for i = 1:size(x, 1)
        p(i,:) = mvncdf(aRegion(1,:), aRegion(2, :), x(i,:), diag(sx)) + mvncdf(bRegion(1,:), bRegion(2, :), x(i,:), diag(sx)) + mvncdf(cRegion(1,:), cRegion(2,:), x(i,:), diag(sx));
%         q(i,:) = mvncdf(dRegion(1,:), dRegion(2, :), x(i,:), diag(sx));
    end 
end
p = 1 - p; 