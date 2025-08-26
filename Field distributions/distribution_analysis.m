opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

tbl = readtable("Weise2025.csv", opts);

SE = tbl.x;
pSE = tbl.y;

clear opts tbl
%%
[SE, ind] = sort(SE);
pSE = pSE(ind);
SE_vec = (0:0.02:5)';
pSE_vec = interp1(SE, pSE, SE_vec, "linear", "extrap");


Q = trapz(SE_vec, pSE_vec);
pSE_vec = pSE_vec / Q;

mean_x = trapz(SE_vec, pSE_vec.*SE_vec);

SEinterp_vec = (SE_vec(1:end-1) + SE_vec(2:end))/2;
pSEinterp_vec = interp1(SE_vec, pSE_vec, SEinterp_vec, "linear", "extrap");

save("pSE_distribution.mat", "SE*", "pSE*");

%%
figure;
plot(SE, pSE, '.k', SE_vec, pSE_vec, 'k-');
hold on;
plot(SEinterp_vec, pSEinterp_vec, 'g--');

axis([0, 5, 0, max(ylim)]);