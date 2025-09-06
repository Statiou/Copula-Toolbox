function make_demo()
% Create a small demo.csv in the examples folder.
rng default
n = 800; R = [1 .6 .3; .6 1 .4; .3 .4 1]; nu = 6;
U = copularnd('t', R, n, nu);
X = [norminv(U(:,1)), norminv(U(:,2)), norminv(U(:,3))];
T = array2table(X, 'VariableNames', {'X1','X2','X3'});
writetable(T, fullfile(fileparts(mfilename('fullpath')), 'demo.csv'));
end
