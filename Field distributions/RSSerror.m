function RSS = RSSerror(params, N, lambda, phi, z_vec, SE_vec, pSEint_vec)

A = sort(abs(params));
% A = abs(params);
E_z = ones(size(z_vec));        % mV/cm
    
for ii = 1 : N
    E_z = E_z + A(ii) * sin( 2*pi * z_vec/lambda(ii) - phi(ii));
end

pSE_actual = histcounts(abs(E_z), SE_vec, 'Normalization', 'pdf')';

RSS = ( sum( (pSEint_vec - pSE_actual).^2, 'all', 'omitnan') );

end