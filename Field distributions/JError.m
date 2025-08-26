function J = JError(lambda, x1, x2, W)   % mean square of absolute error
% x(1) = tau, x(2) = V_max

if nargin < 4
    W = 1;
end

J = sum((x1 * lambda - x2).^2 .* W, 'all', 'omitnan');
if (lambda <= 0)
    J = inf;
end

end
