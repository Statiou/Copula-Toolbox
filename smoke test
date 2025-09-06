function ok = smoke_tests()
% Minimal headless checks for JOSS reviewers (no UI needed).
% MATLAB R2023b + Statistics and Machine Learning Toolbox.

rng default

%% --- Load or synthesize demo data ---
T = [];
if exist('sampledata.m','file') == 2
    T = sampledata();             % should return a table with >=2 cols
    if ~istable(T) || width(T) < 2
        error('sampledata() must return a table with at least 2 columns.');
    end
else
    % fallback: synthesize 2D t-copula and map via marginals
    n  = 600;
    R  = [1 .6; .6 1];
    nu = 6;
    U  = copularnd('t', R, n, nu);
    U  = min(max(U,1e-12),1-1e-12);

    % Map *both* margins from the same U (ensures identical lengths)
    x1 = norminv(U(:,1));           x1 = x1(:);
    x2 = logninv(U(:,2),0,1);       x2 = x2(:);

    X  = [x1, x2];
    T  = array2table(X, 'VariableNames', {'X1','X2'});
end

X = T{:,1:2};

%% --- PIT under parametric marginals ---
pd1 = fitdist(X(:,1),'Normal');
pd2 = fitdist(X(:,2),'Lognormal');
U   = [cdf(pd1,X(:,1)), cdf(pd2,X(:,2))];
U   = min(max(U,1e-12),1-1e-12);

%% --- Fit Gaussian and t copulas ---
[Rg]      = copulafit('Gaussian', U);
[Rt,nu_t] = copulafit('t', U);

%% --- Basic assertions ---
assert(all(size(Rg)==[2 2]) && all(eig(Rg) > 0), 'Rg must be PD 2x2');
assert(nu_t > 2, 't DoF must be >2');

%% --- Densities/logL must be finite ---
Lg = sum(log(max(copulapdf('Gaussian',U,Rg), realmin)));
Lt = sum(log(max(copulapdf('t',U,Rt,nu_t),  realmin)));
assert(isfinite(Lg) && isfinite(Lt), 'log-likelihoods not finite');

%% --- Save a quick artifact (proves plotting works) ---
f = figure('Visible','off'); plotmatrix(X); title('smoke: scatter');
saveas(f, fullfile(pwd,'smoke_scatter.png')); close(f);

disp(table(Lg,Lt));
ok = true;
end
