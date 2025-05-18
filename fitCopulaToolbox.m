% fitCopulaToolbox.m

function fitCopulaToolbox
    %% Splash Screen
    splash = uifigure('Name','Starting...','Position',[500 400 400 200], 'Color','white');
    uilabel(splash, 'Text','Copula Toolbox by Hatzopoulos Petros & Statiou D. Anastasios', 'Position',[30 100 340 40], ...
        'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment','center');
    uilabel(splash, 'Text','Loading, please wait...', 'Position',[30 60 340 30], ...
        'FontSize', 12, 'HorizontalAlignment','center');
    pause(2);
    close(splash);

    %% Main GUI
    f = uifigure('Name','Copula Fitting Toolbox by Xatzopoulos Petros & Statiou. D. Anastasios','Position',[100 100 1000 600]);

    % UI Components
    uilabel(f,'Position',[20 530 100 22],'Text','Copula:');
    ddlCopula = uidropdown(f,'Position',[100 530 200 22], 'Items',{'Gaussian','t','Clayton','Frank','Gumbel'});

    uilabel(f,'Position',[20 490 100 22],'Text','Marginals:');
    ddlMarg = uidropdown(f,'Position',[100 490 200 22],'Items',{'Normal','Lognormal','Exponential'});

    btnLoad = uibutton(f,'Text','Load Data','Position',[320 530 100 22],'ButtonPushedFcn', @(btn,event) loadData());
    uilabel(f,'Position',[20 450 100 22],'Text','Columns:');
    lstCols = uilistbox(f,'Position',[100 310 200 140],'Multiselect','on');

    btnFit = uibutton(f,'Text','Fit Copula','Position',[320 490 110 24],'ButtonPushedFcn', @(btn,event) fitMultivariate());
    btnMatrix = uibutton(f,'Text','Scatter Matrix','Position',[440 490 130 24],'ButtonPushedFcn', @(btn,event) showScatterMatrix());
    btnCompareR = uibutton(f,'Text','Compare R-matrices','Position',[580 490 150 24],'ButtonPushedFcn', @(btn,event) compareCorrelations());
    btnMarginals = uibutton(f,'Text','Marginal Histograms','Position',[740 490 150 24],'ButtonPushedFcn', @(btn,event) showMarginals());
    btnSimulate = uibutton(f,'Text','Simulate from Copula','Position',[900 490 160 24],'ButtonPushedFcn', @(btn,event) simulateFromCopula());
    btnJointPDF = uibutton(f,'Text','Plot Joint PDF (2D)','Position',[320 450 160 24],'ButtonPushedFcn', @(btn,event) plotJointPDF());
    btnAutoSelect = uibutton(f,'Text','Auto-Select Copula','Position',[500 450 160 24],'ButtonPushedFcn', @(btn,event) autoSelectCopula());

    ax = uiaxes(f,'Position',[500 60 460 400]);
    title(ax, 'Correlation Heatmap');

    txtOut = uitextarea(f,'Position',[20 100 460 200],'Editable','off');
    uilabel(f, 'Position', [760 560 200 22], 'Text', 'Statiou D. Anastasios', 'FontAngle', 'italic', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

    userData = struct('Table', [], 'U', [], 'R', [], 'nu', []);

    %% Load Data
    function loadData()
        [file, path] = uigetfile({'*.csv;*.mat','Data Files (*.csv, *.mat)'},'Select Data File');
        if isequal(file,0), return; end
        fullpath = fullfile(path,file);
        [~,~,ext] = fileparts(fullpath);
        switch lower(ext)
            case '.csv'; T = readtable(fullpath);
            case '.mat'; vars = load(fullpath); fn = fieldnames(vars); T = vars.(fn{1});
            otherwise; uialert(f,'Unsupported file format','Error'); return;
        end
        if ~istable(T), T = array2table(T); end
        userData.Table = T;
        lstCols.Items = T.Properties.VariableNames;
    end

    %% Fit Multivariate Copula
    function fitMultivariate()
        cols = lstCols.Value;
        if isempty(cols) || numel(cols) < 2
            uialert(f,'Select at least 2 columns','Error'); return;
        end
        X = userData.Table{:,cols};
        marg = lower(ddlMarg.Value);
        U = zeros(size(X));
        for i = 1:size(X,2)
            pd = fitdist(X(:,i), marg);
            U(:,i) = cdf(pd, X(:,i));
        end

        copulaType = lower(ddlCopula.Value);
        switch copulaType
            case 'gaussian'
                R = copulafit('Gaussian', U);
                nu = [];
            case 't'
                [R, nu] = copulafit('t', U);
                userData.nu = nu;
            case {'clayton','frank','gumbel'}
                if size(U,2) ~= 2
                    uialert(f, 'Clayton, Frank and Gumbel copulas only support 2D data.', 'Copula Restriction');
                    return;
                end
                theta = copulafit(copulaType, U);
                R = []; nu = []; % δεν υπάρχει πίνακας συσχέτισης
                txtOut.Value = sprintf('Fitted bivariate %s copula with parameter θ = %.4f\nMarginals: %s', ...
                       ddlCopula.Value, theta, marg);

                userData.U = U; userData.R = []; userData.theta = theta;
                return;
            otherwise
                uialert(f,'Unsupported copula','Error'); return;
        end

        userData.U = U; userData.R = R;
        txt = sprintf('Fitted multivariate %s copula on %d variables\nMarginals: %s\n', ddlCopula.Value, size(U,2), marg);
        if exist('nu','var') && ~isempty(nu)
            txt = [txt, sprintf('Degrees of freedom (nu): %.4f\n', nu)];
        end
        txt = [txt, evalc('disp(R)')];
        txtOut.Value = txt;
        imagesc(ax, R); colorbar(ax);
        xticks(ax, 1:numel(cols)); yticks(ax, 1:numel(cols));
        xticklabels(ax, cols); yticklabels(ax, cols);
    end

    %% Scatter Matrix
    function showScatterMatrix()
        if isempty(userData.Table) || isempty(lstCols.Value)
            uialert(f,'Load data and select columns first','Error'); return;
        end
        X = userData.Table{:,lstCols.Value};
        figure; plotmatrix(X);
        sgtitle('Scatter Matrix of Selected Variables');
    end

    %% Compare Pearson vs Copula Correlation
    function compareCorrelations()
        if isempty(userData.Table) || isempty(userData.R)
            uialert(f,'Run fitting first','Error'); return;
        end
        X = userData.Table{:,lstCols.Value};
        pearsonR = corr(X);
        copulaR = userData.R;
        figure;
        subplot(1,2,1); imagesc(pearsonR); colorbar; title('Pearson Correlation');
        xticks(1:size(X,2)); yticks(1:size(X,2));
        subplot(1,2,2); imagesc(copulaR); colorbar; title('Copula Implied Correlation');
        xticks(1:size(X,2)); yticks(1:size(X,2));
    end

    %% Show Marginals
    function showMarginals()
        if isempty(userData.Table) || isempty(lstCols.Value)
            uialert(f,'Load data and select columns first','Error'); return;
        end
        cols = lstCols.Value;
        X = userData.Table{:,cols};
        marg = lower(ddlMarg.Value);
        figure;
        for i = 1:size(X,2)
            subplot(ceil(sqrt(size(X,2))), ceil(sqrt(size(X,2))), i);
            histogram(X(:,i),'Normalization','pdf'); hold on;
            pd = fitdist(X(:,i), marg);
            x_range = linspace(min(X(:,i)), max(X(:,i)), 100);
            plot(x_range, pdf(pd, x_range), 'r','LineWidth',1.5);
            title(cols{i});
        end
        sgtitle('Marginal Histograms with Fitted Distributions');
    end

    %% Simulate from Copula
    function simulateFromCopula()
        if isempty(userData.R), uialert(f,'Fit a copula first','Error'); return; end
        n = size(userData.U,1);
        d = size(userData.U,2);
        marg = lower(ddlMarg.Value);
        cols = lstCols.Value;
        if strcmpi(ddlCopula.Value,'t') && isfield(userData,'nu')
            nu = round(userData.nu);
            U_sim = copularnd('t', userData.R, n, nu);
        else
            U_sim = copularnd('Gaussian', userData.R, n);
        end
        X_sim = zeros(size(U_sim));
        for i = 1:d
            x_real = userData.Table{:,cols{i}};
            pd = fitdist(x_real, marg);
            X_sim(:,i) = icdf(pd, U_sim(:,i));
        end
        figure; plotmatrix(X_sim); sgtitle('Simulated Data from Copula');
    end

    %% Plot Joint PDF
    function plotJointPDF()
        cols = lstCols.Value;
        if numel(cols) ~= 2
            uialert(f, 'Select exactly 2 columns for joint PDF plot.', 'Invalid Selection'); return;
        end
        X1 = userData.Table{:,cols{1}};
        X2 = userData.Table{:,cols{2}};
        marg = lower(ddlMarg.Value);
        copulaType = lower(ddlCopula.Value);
        pd1 = fitdist(X1, marg); pd2 = fitdist(X2, marg);
        U1 = cdf(pd1, X1); U2 = cdf(pd2, X2); U = [U1, U2];
        if strcmp(copulaType, 't')
            [R, nu] = copulafit('t', U);
        else
            R = copulafit('Gaussian', U); nu = [];
        end
        x1 = linspace(min(X1), max(X1), 50);
        x2 = linspace(min(X2), max(X2), 50);
        [Xgrid, Ygrid] = meshgrid(x1, x2);
        U1_grid = cdf(pd1, Xgrid); U2_grid = cdf(pd2, Ygrid);
        f1 = pdf(pd1, Xgrid); f2 = pdf(pd2, Ygrid);
        if strcmp(copulaType, 't')
            z1 = tinv(U1_grid(:), nu); z2 = tinv(U2_grid(:), nu);
            Z = [z1 z2];
            c = mvtpdf(Z, R, nu) ./ (tpdf(z1,nu) .* tpdf(z2,nu));
        else
            Z = [norminv(U1_grid(:)), norminv(U2_grid(:))];
            c = mvnpdf(Z, [0 0], R) ./ (normpdf(Z(:,1)) .* normpdf(Z(:,2)));
        end
        c = reshape(c, size(Xgrid));
        jointPDF = c .* f1 .* f2;
        figure; surf(x1, x2, jointPDF, 'EdgeColor','none');
        xlabel(cols{1}); ylabel(cols{2}); zlabel('Joint PDF');
        title(['Joint PDF: ', ddlCopula.Value, ' + ', ddlMarg.Value]);
        view(135, 30); colorbar;
    end

    %% Auto-Select Copula
    function autoSelectCopula()
        cols = lstCols.Value;
        if isempty(cols) || numel(cols) < 2
            uialert(f,'Select at least 2 columns','Error'); return;
        end
        X = userData.Table{:,cols};
        marg = lower(ddlMarg.Value);
        U = zeros(size(X));
        for i = 1:size(X,2)
            pd = fitdist(X(:,i), marg);
            U(:,i) = cdf(pd, X(:,i));
        end

        families = {'Gaussian','t','Clayton','Frank','Gumbel'};
        results = {};

        for k = 1:numel(families)
            try
                fam = families{k};
                if strcmpi(fam, 't')
                    [params, nu] = copulafit(fam, U);
                    logL = sum(log(copulapdf(fam, U, params, nu)));
                    kparams = numel(params) + 1;
                else
                    params = copulafit(fam, U);
                    logL = sum(log(copulapdf(fam, U, params)));
                    kparams = numel(params);
                end
                n = size(U,1);
                AIC = -2*logL + 2*kparams;
                BIC = -2*logL + kparams*log(n);
                results(end+1,:) = {fam, logL, AIC, BIC};
            catch
                results(end+1,:) = {fam, NaN, NaN, NaN};
            end
        end

        T = cell2table(results, 'VariableNames', {'Copula','LogL','AIC','BIC'});

        [~, iAIC] = min(T.AIC);
        [~, iBIC] = min(T.BIC);

       msg = sprintf('Model selection completed:\nBest AIC: %s\nBest BIC: %s', ...
              T.Copula{iAIC}, T.Copula{iBIC});

        uialert(f, msg, 'Copula Selection Results');

        % Πίνακας σε νέο παράθυρο
        figResults = uifigure('Name','Model Selection Table','Position',[300 300 520 200]);
        data = [string(T.Copula), ...
                string(round(T.LogL,3)), ...
                string(round(T.AIC,3)), ...
                string(round(T.BIC,3))];
        uitable(figResults, ...
            'Data', data, ...
            'ColumnName', {'Copula','LogL','AIC','BIC'}, ...
            'Position', [20 20 480 150]);
    end
end
