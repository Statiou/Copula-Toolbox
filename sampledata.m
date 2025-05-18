rng(1); % για αναπαραγωγιμότητα
n = 100;

% Δημιουργία συσχετισμένων μεταβλητών (Gaussian)
mu = [0 0 0 0];
Sigma = [1 0.8 0.5 0.2; 
         0.8 1 0.3 0.1; 
         0.5 0.3 1 0.4; 
         0.2 0.1 0.4 1];

X = mvnrnd(mu, Sigma, n);

% Προσθήκη θετικού μέσου
X = X + 3;

% Δημιουργία πίνακα και αποθήκευση
T = array2table(X, 'VariableNames', {'X1','X2','X3','X4'});
writetable(T, 'sample_data.csv');
disp('Το αρχείο sample_data.csv δημιουργήθηκε.');
