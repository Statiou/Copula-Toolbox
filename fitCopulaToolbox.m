function fitcopulatoolbox
% Copula Toolbox (single-file UI) — extended & fixed
% -------------------------------------------------------------------------
% Features
% - Multivariate Gaussian / t copulas (d>=2)
% - Bivariate Archimedean: Clayton / Frank / Gumbel
% - Marginals: Normal, Lognormal, Exponential, Gamma, tLocationScale,
%              or Empirical Ranks
% - Scatter matrix, marginal histograms, correlation heatmaps
% - Rank heatmaps (Kendall/Spearman)
% - PIT diagnostics
% - Empirical copula contours (2D)
% - Copula PDF/CDF on U (2D)  [PAIR-SPECIFIC REFIT -> robust]
% - Kendall plot (2D)         [PAIR-SPECIFIC REFIT -> robust]
% - Auto-Select copula (AIC/BIC)
% - K-fold cross-validated log-likelihood
% - Simulation (choose nSim) & inverse mapping via chosen marginals
% - 2D joint PDF surface (any copula + chosen marginals)
% - Tail dependence table (λ_L, λ_U)     [SUBSET REFIT FOR ELLIPTICAL]
% - GOF (2D) via Rosenblatt transform + KS tests
% - Rolling Kendall's τ plot
% - Forecast τ(t+1) with GAS(1,1) + split-conformal bands
% - Parametric bootstrap CIs για R/ν/θ
% - 3D scatter: raw (3 cols) ή U-space (2 cols) με Z=c(u) ή C_n(u)
% - Export U/R/results, Save/Load session, Export Axes PNG, Quick Report
%
% Author: Statiou. D. Anastasios.  (MATLAB R2023b)
% -------------------------------------------------------------------------

%% Splash
splash = uifigure('Name','Starting...','Position',[520 420 460 220],'Color','white');
uilabel(splash,'Text','Copula Toolbox (single-file)',...
    'Position',[30 130 400 32],'FontSize',16,'FontWeight','bold','HorizontalAlignment','center');
uilabel(splash,'Text','Statiou. D. Anastasios.',...
    'Position',[30 96 400 28],'FontSize',18,'FontWeight','bold','HorizontalAlignment','center');
uilabel(splash,'Text','Loading, please wait...','Position',[30 64 400 24],...
    'FontSize',12,'HorizontalAlignment','center');
pause(0.9); try, close(splash); end

%% Main GUI
f = uifigure('Name','Copula Toolbox — Statiou. D. Anastasios.', 'Position',[60 50 1320 760]);

% Left panel (controls)
uilabel(f,'Position',[20 710 120 22],'Text','Copula:');
ddlCopula = uidropdown(f,'Position',[120 710 230 24], ...
    'Items',{'Gaussian','t','Clayton (2D)','Frank (2D)','Gumbel (2D)'}, 'Value','Gaussian');

uilabel(f,'Position',[20 680 120 22],'Text','Marginals:');
ddlMarg = uidropdown(f,'Position',[120 680 230 24], ...
    'Items',{'Normal','Lognormal','Exponential','Gamma','tLocationScale','Ranks (empirical)'}, ...
    'Value','Normal');

btnLoad = uibutton(f,'push','Text','Load Data','Position',[370 710 120 26], ...
    'ButtonPushedFcn', @(src,~) loadData(src));
uilabel(f,'Position',[20 650 120 22],'Text','Columns:');
lstCols = uilistbox(f,'Position',[120 510 230 130],'Multiselect','on','Items',{},'Value',string.empty(0,1));

% Numeric indices (optional)
uilabel(f,'Position',[20 480 120 22],'Text','Col idx:');
edtIdx  = uieditfield(f,'text','Position',[120 480 230 24],'Placeholder','e.g. 1,2,3','Value','');

% Extra controls
uilabel(f,'Position',[20 450 120 22],'Text','alpha:');
edtAlpha = uieditfield(f,'numeric','Position',[120 450 80 24],'Value',0.10,'Limits',[0.001 0.5]);
uilabel(f,'Position',[210 450 60 22],'Text','K-fold:');
edtKfold = uieditfield(f,'numeric','Position',[270 450 80 24],'Value',5,'Limits',[2 20],'RoundFractionalValues',true);
uilabel(f,'Position',[20 420 120 22],'Text','Bootstrap B:');
edtBootB = uieditfield(f,'numeric','Position',[120 420 80 24],'Value',200,'Limits',[50 5000],'RoundFractionalValues',true);

% Action buttons row 1
uibutton(f,'push','Text','Fit Copula','Position',[510 710 120 26],...
    'ButtonPushedFcn', @(src,~) fitMultivariate(src));
uibutton(f,'push','Text','Scatter Matrix','Position',[640 710 130 26],...
    'ButtonPushedFcn', @(src,~) showScatterMatrix(src));
uibutton(f,'push','Text','Compare R-matrices','Position',[780 710 150 26],...
    'ButtonPushedFcn', @(src,~) compareCorrelations(src));
uibutton(f,'push','Text','Marginal Histograms','Position',[940 710 150 26],...
    'ButtonPushedFcn', @(src,~) showMarginals(src));
uibutton(f,'push','Text','Auto-Select Copula','Position',[1100 710 150 26],...
    'ButtonPushedFcn', @(src,~) autoSelectCopula(src));

% Action buttons row 2
uibutton(f,'push','Text','Joint PDF (2D)','Position',[510 680 120 26],...
    'ButtonPushedFcn', @(src,~) plotJointPDF(src));
uibutton(f,'push','Text','Tail Dependence','Position',[640 680 130 26],...
    'ButtonPushedFcn', @(src,~) tailDependenceDialog(src));
uibutton(f,'push','Text','GOF (Rosenblatt, 2D)','Position',[780 680 150 26],...
    'ButtonPushedFcn', @(src,~) gofDialog(src));
uibutton(f,'push','Text','Rolling \tau','Position',[940 680 150 26],...
    'ButtonPushedFcn', @(src,~) rollingTauDialog(src));
uibutton(f,'push','Text','Forecast \tau(t+1)','Position',[1100 680 150 26],...
    'ButtonPushedFcn', @(src,~) forecastDialog(src));

% Action buttons row 3
uibutton(f,'push','Text','PIT Diagnostics','Position',[510 650 120 26],...
    'ButtonPushedFcn', @(src,~) pitDiagnostics(src));
uibutton(f,'push','Text','Rank Heatmaps','Position',[640 650 130 26],...
    'ButtonPushedFcn', @(src,~) rankHeatmaps(src));
uibutton(f,'push','Text','Empirical Copula (2D)','Position',[780 650 150 26],...
    'ButtonPushedFcn', @(src,~) empiricalCopulaContours(src));
uibutton(f,'push','Text','Copula PDF on U (2D)','Position',[940 650 150 26],...
    'ButtonPushedFcn', @(src,~) copulaPDFonU(src));
uibutton(f,'push','Text','Copula CDF on U (2D)','Position',[1100 650 150 26],...
    'ButtonPushedFcn', @(src,~) copulaCDFonU(src));

% Action buttons row 4
uibutton(f,'push','Text','Kendall Plot (2D)','Position',[510 620 120 26],...
    'ButtonPushedFcn', @(src,~) kendallPlot2D(src));
uibutton(f,'push','Text','K-fold CV LogLik','Position',[640 620 130 26],...
    'ButtonPushedFcn', @(src,~) kfoldCV(src));
uibutton(f,'push','Text','Param Bootstrap','Position',[780 620 150 26],...
    'ButtonPushedFcn', @(src,~) bootstrapDialog(src));
uibutton(f,'push','Text','3D Scatter (2D/3D)','Position',[940 620 150 26],...
    'ButtonPushedFcn', @(src,~) Uscatter3D(src));
uibutton(f,'push','Text','Conditional Slice (2D)','Position',[1100 620 150 26],...
    'ButtonPushedFcn', @(src,~) conditionalSliceDialog(src));

% Sim/Export row
uilabel(f,'Position',[510 590 40 22],'Text','nSim:');
edtNSim = uieditfield(f,'numeric','Position',[550 590 80 24],'Value',1000,'Limits',[10 Inf],'RoundFractionalValues',true);
uibutton(f,'push','Text','Simulate from Copula','Position',[640 590 170 26],...
    'ButtonPushedFcn', @(src,~) simulateFromCopula(src));
uibutton(f,'push','Text','Export U/R/fit','Position',[820 590 130 26],...
    'ButtonPushedFcn', @(src,~) exportResults(src));
uibutton(f,'push','Text','Save Session','Position',[960 590 120 26],...
    'ButtonPushedFcn', @(src,~) saveSession(src));
uibutton(f,'push','Text','Load Session','Position',[1090 590 120 26],...
    'ButtonPushedFcn', @(src,~) loadSession(src));

% Export graphics / report
uibutton(f,'push','Text','Export Axes PNG','Position',[20 390 160 26],...
    'ButtonPushedFcn', @(src,~) exportAxesPNG(src));
uibutton(f,'push','Text','Quick Report','Position',[190 390 160 26],...
    'ButtonPushedFcn', @(src,~) quickReport(src));

% Right: axes + output
ax = uiaxes(f,'Position',[380 120 900 450]); title(ax, 'Correlation Heatmap');
txtOut = uitextarea(f,'Position',[20 120 340 250],'Editable','off','FontName','Consolas');
uilabel(f,'Position',[980 740 320 22],'Text','Statiou. D. Anastasios.','FontAngle','italic','FontWeight','bold','HorizontalAlignment','right');

% State
userData = struct('Table',[], 'U',[], 'R',[], 'nu',[], 'theta',[], 'family','', 'marg','Normal');

%% ---------- Helpers ----------
    function cols = getSelectedColumns()
        cols = [];
        if ~isempty(strtrim(edtIdx.Value)) && ~isempty(userData.Table)
            txt  = regexprep(edtIdx.Value,'\s+','');
            nums = regexp(txt,'\d+','match');
            if ~isempty(nums)
                idx = unique(str2double(nums));
                idx = idx(~isnan(idx) & idx>=1 & idx<=width(userData.Table));
                if ~isempty(idx)
                    cols = userData.Table.Properties.VariableNames(idx);
                end
            end
        end
        if isempty(cols), cols = lstCols.Value; end
        if isstring(cols) || iscategorical(cols), cols = cellstr(cols); end
        if ischar(cols), cols = {cols}; end
        cols = cols(:);
    end

    function U = makeUfromMarginals(X, margName)
        n = size(X,1); d = size(X,2);
        U = zeros(n,d);
        if strcmpi(margName,'Ranks (empirical)')
            for j=1:d
                [~,r] = sort(X(:,j)); invr = zeros(n,1); invr(r) = 1:n; U(:,j) = invr/(n+1);
            end
        else
            for j=1:d
                pd = fitdist(X(:,j), lower(margName));
                U(:,j) = cdf(pd, X(:,j));
            end
        end
        U = min(max(U, 1e-12), 1-1e-12);
    end

    function assertFitted(src)
        if isempty(userData.U)
            showAlert(src,'Fit a copula first.','Info');
            error('Not fitted.');
        end
    end

%% ---------- Load Data ----------
    function loadData(src)
        try
            [file, path] = uigetfile({'*.csv;*.xlsx;*.xls;*.mat','Data Files (*.csv, *.xlsx, *.xls, *.mat)'}, 'Select Data File');
            if isequal(file,0), return; end
            fullpath = fullfile(path,file); [~,~,ext] = fileparts(fullpath);
            switch lower(ext)
                case '.csv'
                    T = readtable(fullpath,'PreserveVariableNames',true);
                case {'.xlsx','.xls'}
                    T = readtable(fullpath,'PreserveVariableNames',true);
                case '.mat'
                    s = load(fullpath); fn = fieldnames(s); T = s.(fn{1});
                otherwise
                    showAlert(src,'Unsupported file format','Error'); return;
            end
            if ~istable(T), T = array2table(T); end
            names = T.Properties.VariableNames(:);
            if isstring(names) || iscategorical(names), names = cellstr(names); end
            names = cellfun(@char, names, 'UniformOutput', false);
            userData.Table = T; lstCols.Value = string.empty(0,1); lstCols.Items = names;
            if numel(names) >= 2, lstCols.Value = names(1:2); else, lstCols.Value = names(1); end
            edtIdx.Value = '';
            txtOut.Value = ""; userData.U=[]; userData.R=[]; userData.nu=[]; userData.theta=[]; userData.family='';
            cla(ax); title(ax,'Correlation Heatmap');
        catch ME
            showAlert(src, ME.message, 'Load Error');
        end
    end

%% ---------- Fit ----------
    function fitMultivariate(src)
        cols = getSelectedColumns();
        if isempty(cols) || numel(cols) < 2
            showAlert(src,'Select at least 2 columns','Error'); return;
        end
        X = userData.Table{:,cols};
        marg = ddlMarg.Value; userData.marg = marg;
        U = makeUfromMarginals(X, marg);
        copulaType = ddlCopula.Value;
        try
            switch lower(copulaType)
                case 'gaussian'
                    R = copulafit('Gaussian', U);
                    R = nearestPD(R); nu = []; theta = []; famKey = 'Gaussian';
                case 't'
                    [R, nu] = copulafit('t', U);
                    R = nearestPD(R); theta = []; famKey = 't';
                case 'clayton (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Clayton supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Clayton', U); R=[]; nu=[]; famKey='Clayton';
                case 'frank (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Frank supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Frank', U); R=[]; nu=[]; famKey='Frank';
                case 'gumbel (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Gumbel supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Gumbel', U); R=[]; nu=[]; famKey='Gumbel';
                otherwise
                    showAlert(src,'Unsupported copula','Error'); return;
            end
        catch ME
            showAlert(src, ME.message, 'Fit Error'); return;
        end

        userData.U = U; userData.R = R; userData.theta = theta; userData.nu = nu; userData.family = famKey;
        txt = sprintf('Fitted %s copula on %d variables\nMarginals: %s\n', famKey, size(U,2), marg);
        if ~isempty(nu), txt = [txt, sprintf('Degrees of freedom (nu): %.4f\n', nu)]; end
        if ~isempty(R)
            txt = [txt, 'R =', newline, evalc('disp(R)')];
            imagesc(ax, R); colorbar(ax);
            xticks(ax, 1:numel(cols)); yticks(ax, 1:numel(cols)); xticklabels(ax, cols); yticklabels(ax, cols);
            title(ax,'Copula Correlation (R)');
        else
            cla(ax); title(ax, sprintf('%s fitted (2D)', famKey));
        end
        txtOut.Value = txt;
    end

%% ---------- Plots / Analysis ----------
    function showScatterMatrix(src)
        if isempty(userData.Table)
            showAlert(src,'Load data first','Error'); return;
        end
        cols = getSelectedColumns();
        if numel(cols) < 2, showAlert(src,'Select ≥2 columns','Error'); return; end
        X = userData.Table{:,cols};
        figure('Name','Scatter Matrix'); plotmatrix(X); sgtitle('Scatter Matrix of Selected Variables');
    end

    function compareCorrelations(src)
        if isempty(userData.Table)
            showAlert(src,'Load data first','Error'); return;
        end
        if any(strcmpi(userData.family,{'Clayton','Frank','Gumbel'}))
            showAlert(src,'Compare R available only for Gaussian/t.','Info'); return;
        end
        cols = getSelectedColumns(); 
        if numel(cols)<2, showAlert(src,'Select ≥2 columns','Error'); return; end
        X = userData.Table{:,cols};
        U = makeUfromMarginals(X, ddlMarg.Value); % recompute for subset
        try
            switch userData.family
                case 't'
                    [Rsub,~] = copulafit('t', U);
                otherwise
                    Rsub = copulafit('Gaussian', U);
            end
        catch
            Rsub = corr(norminv(U), 'rows','pairwise'); % fallback
        end
        pearsonR = corr(X, 'type','Pearson','rows','pairwise');
        figure('Name','Correlation Comparison');
        subplot(1,2,1); imagesc(pearsonR); colorbar; title('Pearson Correlation');
        xticks(1:size(X,2)); yticks(1:size(X,2)); xticklabels(cols); yticklabels(cols);
        subplot(1,2,2); imagesc(Rsub); colorbar; title('Copula (R) [subset]');
        xticks(1:size(X,2)); yticks(1:size(X,2)); xticklabels(cols); yticklabels(cols);
    end

    function showMarginals(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if isempty(cols), showAlert(src,'Select columns','Error'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value;
        figure('Name','Marginals'); m = size(X,2); r = ceil(sqrt(m)); c = ceil(m/r);
        for i = 1:m
            subplot(r,c,i); histogram(X(:,i),'Normalization','pdf'); hold on;
            if ~strcmpi(marg,'Ranks (empirical)')
                pd = fitdist(X(:,i), lower(marg)); x_range = linspace(min(X(:,i)), max(X(:,i)), 200);
                plot(x_range, pdf(pd, x_range), 'LineWidth',1.4);
            end
            title(cols{i},'Interpreter','none');
        end
        sgtitle(sprintf('Marginal Histograms (marginals: %s)', marg));
    end

    function rankHeatmaps(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2','Error'); return; end
        X = userData.Table{:,cols};
        K = corr(X,'type','Kendall','rows','pairwise');
        S = corr(X,'type','Spearman','rows','pairwise');
        figure('Name','Rank Correlations');
        subplot(1,2,1); imagesc(K); colorbar; title('Kendall \tau');
        xticks(1:numel(cols)); yticks(1:numel(cols)); xticklabels(cols); yticklabels(cols);
        subplot(1,2,2); imagesc(S); colorbar; title('Spearman \rho');
        xticks(1:numel(cols)); yticks(1:numel(cols)); xticklabels(cols); yticklabels(cols);
    end

    function pitDiagnostics(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)<1, showAlert(src,'Select columns','Error'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
        m = size(U,2); r = ceil(sqrt(m)); c = ceil(m/r);
        figure('Name','PIT Diagnostics (U~U[0,1])');
        for j=1:m
            subplot(r,c,j);
            histogram(U(:,j), 'Normalization','pdf'); hold on;
            x = linspace(0,1,200); plot(x, ones(size(x)), 'LineWidth',1.2);
            title(sprintf('%s (PIT)', cols{j}), 'Interpreter','none'); xlim([0 1]);
        end
        sgtitle('PIT histograms: good marginals => approx. flat at 1');
    end

    function empiricalCopulaContours(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Empirical Copula'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value); U1=U(:,1); U2=U(:,2);
        g = linspace(0,1,60); [G1,G2] = meshgrid(g,g);
        Cemp = zeros(size(G1));
        for i=1:numel(G1)
            Cemp(i) = mean(U1<=G1(i) & U2<=G2(i));
        end
        figure('Name','Empirical Copula (CDF)');
        contourf(g,g,Cemp, 12, 'LineStyle','none'); colorbar;
        xlabel('u'); ylabel('v'); title(sprintf('Empirical Copula C_n(u,v): %s-%s', cols{1}, cols{2}));
    end

    % -------- FIXED: pair-specific refit, no logical indexing into userData.U
    function copulaPDFonU(src)
        assertFitted(src);
        cols = getSelectedColumns();
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','PDF on U'); return; end
        Xpair = userData.Table{:,cols};
        Upair = makeUfromMarginals(Xpair, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', Upair); param = {R2, max(2,round(nu2))};
            case 'Gaussian'
                R2 = copulafit('Gaussian', Upair); param = {R2};
            otherwise
                th2 = copulafit(fam, Upair); param = {th2};
        end
        g = linspace(0,1,80); [G1,G2] = meshgrid(g,g);
        switch fam
            case 't',        Cpdf = copulapdf('t', [G1(:) G2(:)], param{1}, param{2});
            case 'Gaussian', Cpdf = copulapdf('Gaussian', [G1(:) G2(:)], param{1});
            otherwise,       Cpdf = copulapdf(fam, [G1(:) G2(:)], param{1});
        end
        Cpdf = reshape(Cpdf,size(G1));
        figure('Name','Copula PDF on U'); surf(g,g,Cpdf,'EdgeColor','none'); view(135,30); colorbar;
        xlabel('u'); ylabel('v'); zlabel('c(u,v)');
        title(sprintf('%s copula PDF on U — pair: %s-%s', fam, cols{1}, cols{2}));
    end

    % -------- also pair-specific refit for CDF on U
    function copulaCDFonU(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','CDF on U'); return; end
        Xpair = userData.Table{:,cols}; Upair = makeUfromMarginals(Xpair, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', Upair); param = {R2, max(2,round(nu2))};
            case 'Gaussian'
                R2 = copulafit('Gaussian', Upair); param = {R2};
            otherwise
                th2 = copulafit(fam, Upair); param = {th2};
        end
        g = linspace(0,1,80); [G1,G2] = meshgrid(g,g);
        switch fam
            case 't',        Ccdf = copulacdf('t', [G1(:) G2(:)], param{1}, param{2});
            case 'Gaussian', Ccdf = copulacdf('Gaussian', [G1(:) G2(:)], param{1});
            otherwise,       Ccdf = copulacdf(fam, [G1(:) G2(:)], param{1});
        end
        Ccdf = reshape(Ccdf,size(G1));
        figure('Name','Copula CDF on U'); contourf(g,g,Ccdf, 12, 'LineStyle','none'); colorbar;
        xlabel('u'); ylabel('v'); title(sprintf('%s copula CDF on U — pair: %s-%s', fam, cols{1}, cols{2}));
    end

    % -------- pair-specific Kendall plot
    function kendallPlot2D(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Kendall Plot'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', U); W = copulacdf('t',U, R2, max(2,round(nu2)));
            case 'Gaussian'
                R2 = copulafit('Gaussian', U); W = copulacdf('Gaussian',U, R2);
            otherwise
                th2 = copulafit(fam, U); W = copulacdf(fam, U, th2);
        end
        [femp,xx] = ecdf(W);
        t = linspace(0,1,200); K0 = t - t.*log(max(t,eps));
        figure('Name','Kendall Plot');
        plot(xx,femp,'LineWidth',1.6); hold on; plot(t,K0,'--','LineWidth',1.2); grid on;
        xlabel('t'); ylabel('K(t)'); legend('Empirical K','Independence K_0(t)=t-t\log t','Location','best');
        title(sprintf('Kendall plot: %s-%s (%s)', cols{1}, cols{2}, fam));
    end

%% ---------- Tail dependence ----------
    function tailDependenceDialog(src)
        cols = getSelectedColumns(); if numel(cols) < 2, showAlert(src,'Select ≥2 columns','Tail'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','Tail'); return; end

        fam = userData.family;
        if any(strcmpi(fam,{'Clayton','Frank','Gumbel'})) && numel(cols)~=2
            showAlert(src,'Tail dep. for Archimedeans is 2D-only here. Select 2 columns.','Tail'); return;
        end

        % For elliptical, recompute R on the selected subset (robust)
        if any(strcmpi(fam,{'Gaussian','t'}))
            X = userData.Table{:,cols};
            U = makeUfromMarginals(X, ddlMarg.Value);
            try
                if strcmpi(fam,'t')
                    [Rsub,nu2] = copulafit('t', U);
                else
                    Rsub = copulafit('Gaussian', U);
                end
            catch
                Rsub = corr(norminv(U),'rows','pairwise');
            end
            names = cols; d = numel(cols); pairs = nchoosek(1:d,2);
            out = cell(size(pairs,1), 4);
            for i=1:size(pairs,1)
                a = pairs(i,1); b = pairs(i,2); rho = Rsub(a,b);
                if strcmpi(fam,'Gaussian')
                    lamL=0; lamU=0;
                else
                    nuv = max(2,round(nu2));
                    tval = -sqrt((nuv+1)*(1-rho)/(1+rho));
                    lam = 2*tcdf(tval, nuv+1);
                    lamL = lam; lamU = lam;
                end
                out(i,:) = {sprintf('%s-%s',names{a},names{b}), rho, lamL, lamU};
            end
        else
            theta = userData.theta; switch fam
                case 'Clayton', lamL = 2^(-1/theta); lamU = 0;
                case 'Gumbel',  lamL = 0; lamU = 2 - 2^(1/theta);
                otherwise,      lamL = 0; lamU = 0; % Frank
            end
            out = {sprintf('%s-%s',cols{1},cols{2}), NaN, lamL, lamU};
        end
        fig = uifigure('Name','Tail Dependence','Position',[360 360 560 260]);
        uitable(fig,'Data',out,'ColumnName',{'Pair','rho','lambda_L','lambda_U'},'Position',[20 20 520 210]);
    end

%% ---------- GOF (2D Rosenblatt) ----------
    function gofDialog(src)
        cols = getSelectedColumns();
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns.','GOF'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','GOF'); return; end
        runGOF(src, cols);
    end

    function runGOF(src, cols)
        X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        fam = userData.family; epsFD = 1e-4;
        switch fam
            case 'Gaussian'
                R = copulafit('Gaussian', U);
            case 't'
                [R,nu] = copulafit('t', U); %#ok<NASGU>
            otherwise
                theta = copulafit(fam, U); %#ok<NASGU>
        end
        V1 = U(:,1); V2 = nan(size(V1));
        for i=1:size(U,1)
            u1 = U(i,1); u2 = U(i,2);
            switch fam
                case 'Gaussian'
                    Cplus = copulacdf('Gaussian',[min(u1+epsFD,1) u2], R);
                    Cminus= copulacdf('Gaussian',[max(u1-epsFD,0) u2], R);
                case 't'
                    Cplus = copulacdf('t',[min(u1+epsFD,1) u2], R, nu);
                    Cminus= copulacdf('t',[max(u1-epsFD,0) u2], R, nu);
                otherwise
                    Cplus = copulacdf(fam,[min(u1+epsFD,1) u2], theta);
                    Cminus= copulacdf(fam,[max(u1-epsFD,0) u2], theta);
            end
            V2(i) = max(0,min(1,(Cplus - Cminus)/(2*epsFD)));
        end
        [~,p1] = kstest(V1,[sort(V1) ( (1:numel(V1))'/numel(V1) )]);
        [~,p2] = kstest(V2,[sort(V2) ( (1:numel(V2))'/numel(V2) )]);
        pR = 1 - abs(corr(V1,V2,'type','Spearman','rows','complete'));
        msg = sprintf('GOF (Rosenblatt 2D)\nKS p(V1~U): %.3f\nKS p(V2~U): %.3f\nIndependence proxy p: %.3f (higher ~ better)', p1,p2,pR);
        showAlert(src, msg, 'GOF results');
    end

%% ---------- Rolling Kendall tau ----------
    function rollingTauDialog(src)
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns.','Rolling'); return; end
        dlg = uifigure('Name','Rolling \tau settings','Position',[260 260 340 160]);
        uilabel(dlg,'Text','Window W:','Position',[20 90 100 22],'HorizontalAlignment','right');
        edW = uieditfield(dlg,'numeric','Position',[130 90 160 22],'Value',60,'Limits',[10 Inf],'RoundFractionalValues',true);
        uibutton(dlg,'push','Text','Run','Position',[110 40 120 26], 'ButtonPushedFcn', @(~,~) runRollingTau(dlg, round(edW.Value)));
    end

    function runRollingTau(dlg, W)
        try
            cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg); n=size(U,1);
            tau = nan(n,1); for t=W:n, tau(t)=corr(U(t-W+1:t,1),U(t-W+1:t,2),'type','Kendall','rows','complete'); end
            figure('Name',sprintf('Rolling \tau (W=%d)',W),'Color','w'); plot(tau,'-','LineWidth',1.3); grid on; ylim([-1 1]); yline(0,'k:');
            xlabel('t'); ylabel('Kendall''s \tau'); title(sprintf('Rolling \tau for %s-%s (W=%d)', cols{1}, cols{2}, W));
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Rolling Error');
        end
    end

%% ---------- Forecast UI (GAS) ----------
    function forecastDialog(src)
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns for forecasting (e.g., X1,X2).','Forecast'); return; end
        dlg = uifigure('Name','Forecast \tau(t+1) settings','Position',[240 240 380 240]);
        uilabel(dlg,'Text','Family:','Position',[20 180 120 22],'HorizontalAlignment','right');
        ddFam = uidropdown(dlg,'Position',[150 180 200 22],'Items',{'Gaussian','t'},'Value','Gaussian');
        uilabel(dlg,'Text','W (rolling \tau):','Position',[20 140 120 22],'HorizontalAlignment','right');
        edW = uieditfield(dlg,'numeric','Position',[150 140 200 22],'Value',60,'Limits',[10 Inf],'RoundFractionalValues',true);
        uilabel(dlg,'Text','\alpha (1-coverage):','Position',[20 100 120 22],'HorizontalAlignment','right');
        edA = uieditfield(dlg,'numeric','Position',[150 100 200 22],'Value',max(0.01,min(0.30,edtAlpha.Value)),'Limits',[0.01 0.3]);
        uilabel(dlg,'Text','Calibration split:','Position',[20 60 120 22],'HorizontalAlignment','right');
        edC = uieditfield(dlg,'numeric','Position',[150 60 200 22],'Value',0.20,'Limits',[0.05 0.9]);
        uibutton(dlg,'push','Text','Run Forecast','Position',[120 20 140 26], ...
            'ButtonPushedFcn', @(~,~)runForecast(dlg, ddFam.Value, round(edW.Value), edA.Value, edC.Value));
    end

    function runForecast(dlg, fam, W, alpha, calFrac)
        try
            cols = getSelectedColumns(); if numel(cols)~=2, showAlert(dlg,'Need exactly 2 columns.','Forecast'); return; end
            X = userData.Table{:,cols}; marg = ddlMarg.Value; out = forecastTauGAS_core(X, marg, W, fam, alpha, calFrac);
            n = size(X,1); tt = (1:n)'; m = ~isnan(out.tau_pred) & ~isnan(out.tau_real);
            fig = figure('Name',sprintf('\tau Forecast (%s, W=%d, \alpha=%.2f)', fam, W, alpha),'Color','w'); %#ok<NASGU>
            plot(tt(m), out.tau_real(m), '-', 'LineWidth',1.2); hold on;
            plot(tt(m), out.tau_pred(m), '--', 'LineWidth',1.6);
            plot(tt(m), out.loA(m), '-.', 'LineWidth',1.0); plot(tt(m), out.hiA(m), '-.', 'LineWidth',1.0);
            plot(tt(m), out.lo90(m), ':', 'LineWidth',1.0); plot(tt(m), out.hi90(m), ':', 'LineWidth',1.0);
            plot(tt(m), out.lo95(m), ':', 'LineWidth',0.9); plot(tt(m), out.hi95(m), ':', 'LineWidth',0.9);
            yline(0,'k:'); grid on; xlabel('t'); ylabel('Kendall''s \tau');
            legend('\tau realized','\tau^ (GAS)','(1-\alpha) lo','(1-\alpha) hi','90% lo','90% hi','95% lo','95% hi','Location','best');
            title(sprintf('%s + GAS | W=%d | \alpha=%.2f | cal=%.2f', fam, W, alpha, calFrac)); hold off;
            L = min(50, sum(m)); idx = find(m, L, 'last');
            fig2 = uifigure('Name','Forecast table','Position',[300 300 720 300]); %#ok<NASGU>
            uitable(fig2,'Data', [num2cell(idx(:)), ...
                                  num2cell(out.tau_real(idx)), num2cell(out.tau_pred(idx)), ...
                                  num2cell(out.loA(idx)), num2cell(out.hiA(idx)), ...
                                  num2cell(out.lo90(idx)), num2cell(out.hi90(idx)), ...
                                  num2cell(out.lo95(idx)), num2cell(out.hi95(idx))], ...
                           'ColumnName', {'t','tau_real','tau_pred','loA','hiA','lo90','hi90','lo95','hi95'}, ...
                           'Position',[20 20 680 250]);
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Forecast Error');
        end
    end

%% ---------- Forecast core (nested) ----------
    function out = forecastTauGAS_core(X, marg, W, family, alpha, calFrac)
        n = size(X,1);
        if n < W+150
            error('Need at least ~%d observations (got %d).', W+150, n);
        end
        U = makeUfromMarginals(X, marg);
        Z = [norminv(U(:,1)), norminv(U(:,2))]; Z(~isfinite(Z))=0; Xex = abs(Z);
        tau_real = nan(n,1);
        for t=W:n, tau_real(t) = corr(U(t-W+1:t,1), U(t-W+1:t,2), 'type','Kendall','rows','complete'); end
        nu = [];
        if strcmpi(family,'t')
            try
                [~, nuHat] = copulafit('t', U(max(1,W):end,:)); nu = max(3, min(100, round(nuHat)));
            catch, nu = 8; end
        end
        idxStart = W+1; idxEnd = n-1; Ttot = idxEnd-idxStart+1; Tcal = max(50, round(calFrac*Ttot)); Ttrain = Ttot - Tcal; trIdx = idxStart:(idxStart+Ttrain-1);
        muX = mean(Xex(trIdx,:),1,'omitnan'); sX  = std(Xex(trIdx,:),0,1,'omitnan'); sX(sX==0)=1; XexS = (Xex - muX)./sX;
        tau0 = tau_real(idxStart-1);
        if isnan(tau0), tau0 = corr(U(max(1,idxStart-W):idxStart-1,1), U(max(1,idxStart-W):idxStart-1,2), 'type','Kendall','rows','complete'); end
        tau0 = clamp(tau0,-0.999,0.999); rho0 = sin(pi*tau0/2); f0 = atanh(clamp(rho0,-0.9999,0.9999));
        theta0 = [0, log(0.1), log(0.9/0.1), 0, 0];
        obj = @(th) gas_negloglik(th, U, XexS, trIdx, f0, family, nu);
        opts = optimset('Display','off','MaxFunEvals',5e4,'MaxIter',5e4);
        thetaHat = fminsearch(obj, theta0, opts); %#ok<NASGU>
        [~, outTrain] = gas_negloglik(thetaHat, U, XexS, trIdx, f0, family, nu); params = outTrain.params;
        [tau_pred, ~] = predict_tau(params, U, XexS, tau_real, idxStart, idxEnd, family, nu);
        mask = ~isnan(tau_pred) & ~isnan(tau_real); idxAll = find(mask & ( (1:n)'>=idxStart & (1:n)'<=idxEnd ));
        Tcal = max(50, round(calFrac*numel(idxAll))); calIdx = idxAll(end-Tcal+1:end);
        resid = abs(tau_real(calIdx) - tau_pred(calIdx)); resid = resid(isfinite(resid));
        if isempty(resid), qA=0.08; q90=0.10; q95=0.14; else, qA=quantile(resid,1-alpha); q90=quantile(resid,0.90); q95=quantile(resid,0.95); end
        loA  = clamp(tau_pred - qA,  -1, 1);  hiA  = clamp(tau_pred + qA,  -1, 1);
        lo90 = clamp(tau_pred - q90, -1, 1);  hi90 = clamp(tau_pred + q90, -1, 1);
        lo95 = clamp(tau_pred - q95, -1, 1);  hi95 = clamp(tau_pred + q95, -1, 1);
        out = struct('tau_real',tau_real,'tau_pred',tau_pred,'loA',loA,'hiA',hiA,'lo90',lo90,'hi90',hi90,'lo95',lo95,'hi95',hi95,'params',params,'nu',nu,'trainIdx',trIdx);
    end

%% ---------- GAS internals ----------
    function [NLL, out] = gas_negloglik(thetaU, U, Xex, trIdx, f0, fam, nu)
        % thetaU = [w, log(a), logit(b), g1, g2]
        w  = thetaU(1); a  = exp(thetaU(2)); b  = 1/(1+exp(-thetaU(3))); g1 = thetaU(4); g2 = thetaU(5); gam = [g1; g2];
        epsFD = 1e-4; n = size(U,1); f = nan(n,1); rho = nan(n,1); NLL = 0; f(trIdx(1)-1) = f0; rho(trIdx(1)-1) = tanh(f0);
        for t = trIdx(1):trIdx(end)
            xlag  = (t-1>=1) * Xex(max(t-1,1),:).'; fpred = w + b*f(t-1) + gam.'*xlag; rho_t = tanh(fpred);
            ll_t  = logCopulaPdf(U(t,:), rho_t, fam, nu); dldr  = dlogc_drho_fd(U(t,:), rho_t, fam, nu, epsFD); s_t = dldr * sech2(fpred);
            f(t)  = fpred + a * s_t; rho(t)= tanh(f(t)); NLL   = NLL - ll_t;
        end
        out = struct('params',struct('w',w,'a',a,'b',b,'g',gam),'f',f,'rho',rho);
    end

    function [tau_pred, fSeries] = predict_tau(params, U, Xex, tau_real, idxStart, idxEnd, fam, nu)
        n = size(U,1); w=params.w; a=params.a; b=params.b; gam=params.g; f = nan(n,1);
        tau0 = tau_real(idxStart-1);
        if isnan(tau0)
            tau0 = corr(U(max(1,idxStart-60):idxStart-1,1), U(max(1,idxStart-60):idxStart-1,2),'type','Kendall','rows','complete');
        end
        tau0 = clamp(tau0,-0.999,0.999); rho0 = sin(pi*tau0/2); f0 = atanh(clamp(rho0,-0.9999,0.9999)); f(idxStart-1)=f0;
        for t = idxStart:idxEnd
            xlag  = (t-1>=1) * Xex(max(t-1,1),:).'; f_pred= w + b*f(t-1) + gam.' * xlag; rho_t = tanh(f_pred);
            dldr  = dlogc_drho_fd(U(t,:), rho_t, fam, nu, 1e-4); s_t   = dldr * sech2(f_pred);
            f(t)  = f_pred + a * s_t;
        end
        tau_pred = nan(n,1);
        for t = idxStart:idxEnd
            xnow = Xex(t,:).'; f_for = w + b * f(t) + gam.' * xnow; rho_for = tanh(f_for);
            tau_pred(t+1) = (2/pi)*asin(clamp(rho_for,-0.9999,0.9999));
        end
        fSeries=f;
    end

    function val = logCopulaPdf(urow, rho, fam, nu)
        u = max(min(urow(:)', 1-1e-10), 1e-10); R = [1 rho; rho 1];
        if strcmpi(fam,'t')
            if isempty(nu), nu = 8; end
            c = copulapdf('t', u, R, nu);
        else
            c = copulapdf('Gaussian', u, R);
        end
        val = log(max(c, realmin));
    end

    function d = dlogc_drho_fd(urow, rho, fam, nu, epsFD)
        r1 = clamp(rho+epsFD, -0.9999, 0.9999); r2 = clamp(rho-epsFD, -0.9999, 0.9999);
        l1 = logCopulaPdf(urow, r1, fam, nu);  l2 = logCopulaPdf(urow, r2, fam, nu); d  = (l1 - l2) / (2*epsFD);
    end

    function y = sech2(x), y = 1./cosh(x).^2; end
    function z = clamp(x,a,b), z = max(a, min(b, x)); end

%% ---------- Advanced: K-fold CV ----------
    function kfoldCV(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2','CV'); return; end
        X = userData.Table{:,cols};
        U = makeUfromMarginals(X, ddlMarg.Value);  % always recompute for subset
        K = max(2, round(edtKfold.Value)); fams = {'Gaussian','t'};
        if size(U,2)==2, fams = [fams, {'Clayton','Frank','Gumbel'}]; end
        rng('default'); idx = crossvalind('Kfold', size(U,1), K);
        meanLL = nan(1,numel(fams));
        for fi=1:numel(fams)
            fam = fams{fi}; ll = zeros(K,1);
            for k=1:K
                tr = idx~=k; te = idx==k; Utr = U(tr,:); Ute = U(te,:);
                try
                    switch fam
                        case 'Gaussian'
                            R = copulafit('Gaussian', Utr);
                            ll(k) = sum(log(max(copulapdf('Gaussian', Ute, R), realmin)));
                        case 't'
                            [R,nu] = copulafit('t', Utr);
                            ll(k) = sum(log(max(copulapdf('t', Ute, R, max(2,round(nu))), realmin)));
                        otherwise
                            theta = copulafit(fam, Utr);
                            ll(k) = sum(log(max(copulapdf(fam, Ute, theta), realmin)));
                    end
                catch
                    ll(k) = -Inf;
                end
            end
            meanLL(fi) = mean(ll);
        end
        figure('Name','K-fold CV Log-Likelihood'); bar(meanLL);
        set(gca,'XTickLabel',fams,'XTick',1:numel(fams)); ylabel('Mean held-out loglik'); grid on;
        title(sprintf('K-fold (K=%d) CV log-likelihood', K));
    end

%% ---------- Advanced: Bootstrap ----------
    function bootstrapDialog(src)
        assertFitted(src);
        B = max(50, round(edtBootB.Value));
        choice = questdlg(sprintf('Run parametric bootstrap with B=%d?',B), 'Bootstrap','Run','Cancel','Run');
        if ~strcmpi(choice,'Run'), return; end
        bootstrapRun(src, B);
    end

    function bootstrapRun(src, B)
        assertFitted(src); fam = userData.family; U = userData.U; d=size(U,2);
        switch fam
            case {'Gaussian','t'}
                Rhat = userData.R; nu = userData.nu; if isempty(Rhat), showAlert(src,'No R','Bootstrap'); return; end
                Rhats = zeros([size(Rhat) B]); nuv = nan(B,1);
                for b=1:B
                    try
                        if strcmpi(fam,'Gaussian')
                            Ub = copularnd('Gaussian', Rhat, size(U,1));
                            Rb = copulafit('Gaussian', Ub);
                            Rhats(:,:,b) = Rb;
                        else
                            Ub = copularnd('t', Rhat, size(U,1), max(2,round(nu)));
                            [Rb, nub] = copulafit('t', Ub);
                            Rhats(:,:,b) = Rb; nuv(b)=nub;
                        end
                    catch
                        Rhats(:,:,b) = nan(size(Rhat)); nuv(b)=nan;
                    end
                end
                pairs = nchoosek(1:d,2); out = cell(size(pairs,1),5);
                for i=1:size(pairs,1)
                    a=pairs(i,1); b=pairs(i,2);
                    v = squeeze(Rhats(a,b,:)); v=v(isfinite(v));
                    if isempty(v), lo=NaN; hi=NaN; med=NaN; else
                        pr = prctile(v,[2.5 50 97.5]); lo=pr(1); med=pr(2); hi=pr(3);
                    end
                    out(i,:)={sprintf('R(%d,%d)',a,b), userData.R(a,b), lo, med, hi};
                end
                fig = uifigure('Name','Bootstrap CI for R','Position',[300 300 640 280]);
                uitable(fig,'Data',out,'ColumnName',{'Entry','Rhat','lo2.5%','median','hi97.5%'},'Position',[20 20 600 240]);
                if strcmpi(fam,'t')
                    nv = nuv(isfinite(nuv)); if ~isempty(nv)
                        pr = prctile(nv,[2.5 50 97.5]);
                        showAlert(src,sprintf('nu: %.2f (2.5%%=%.2f, 97.5%%=%.2f)',userData.nu,pr(1),pr(3)),'Bootstrap nu');
                    end
                end
            otherwise
                theta = userData.theta; if isempty(theta), showAlert(src,'No theta','Bootstrap'); return; end
                TH = nan(B,1);
                for b=1:B
                    try
                        Ub = copularnd(fam, theta, size(U,1));
                        thb = copulafit(fam, Ub); TH(b)=thb;
                    catch
                        TH(b)=nan;
                    end
                end
                v = TH(isfinite(TH)); if isempty(v), pr=[NaN NaN NaN]; else, pr = prctile(v,[2.5 50 97.5]); end
                showAlert(src,sprintf('theta: %.4f (2.5%%=%.4f, 97.5%%=%.4f)',theta,pr(1),pr(3)),'Bootstrap theta');
        end
    end

%% ---------- Advanced: 3D scatter ----------
    function Uscatter3D(src)
        if isempty(userData.Table), showAlert(src,'Load data first','3D'); return; end
        cols = getSelectedColumns();
        if numel(cols) >= 3
            X = userData.Table{:,cols(1:3)};
            figure('Name','3D Scatter'); scatter3(X(:,1),X(:,2),X(:,3),20,X(:,3),'filled'); grid on;
            xlabel(cols{1}); ylabel(cols{2}); zlabel(cols{3}); title('3D Scatter (raw)'); colorbar; return;
        elseif numel(cols) == 2
            try
                fig = ancestor(src,'figure'); if isempty(fig) || ~isvalid(fig), fig = f; end
                choice = uiconfirm(fig,'Use Z = copula density c(u) (needs fit), or Z = empirical copula C_n(u)?', ...
                                   'Choose Z', 'Options',{'c(u) (fitted)','C_n(u)'}, ...
                                   'DefaultOption',1,'CancelOption',2);
            catch
                choice = 'C_n(u)';
            end
            X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
            U1 = U(:,1); U2 = U(:,2); Z = [];
            if strcmp(choice,'c(u) (fitted)')
                if isempty(userData.family)
                    showAlert(src,'Fit a copula first, then try again (for c(u)).','3D'); return;
                end
                fam = userData.family;
                try
                    switch fam
                        case 'Gaussian'
                            R = copulafit('Gaussian',[U1 U2]); Z = copulapdf('Gaussian',[U1 U2],R);
                        case 't'
                            [R,nu] = copulafit('t',[U1 U2]); Z = copulapdf('t',[U1 U2],R,max(2,round(nu)));
                        otherwise
                            th = copulafit(fam,[U1 U2]); Z = copulapdf(fam,[U1 U2],th);
                    end
                catch
                    showAlert(src,'Could not compute copula density. Falling back to C_n(u).','3D');
                    choice = 'C_n(u)';
                end
            end
            if isempty(Z) || strcmp(choice,'C_n(u)')
                n = numel(U1); Z = zeros(n,1);
                for i=1:n, Z(i) = mean( (U1<=U1(i)) & (U2<=U2(i)) ); end
            end
            figure('Name','U-Scatter 3D (derived Z)'); scatter3(U1,U2,Z,18,Z,'filled'); grid on; xlim([0 1]); ylim([0 1]);
            xlabel('U1'); ylabel('U2'); if strcmp(choice,'C_n(u)'), zlabel('Z = C_n(u)'); else, zlabel('Z = c(u)'); end
            title(sprintf('3D on U for %s-%s', cols{1}, cols{2})); colorbar; return;
        else
            showAlert(src,'Select 2 or 3 columns.','3D'); return;
        end
    end

%% ---------- Advanced: Conditional slice ----------
    function conditionalSliceDialog(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Conditional'); return; end
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Conditional'); return; end
        dlg = uifigure('Name','Conditional slice U2|U1=u','Position',[260 260 380 180]);
        uilabel(dlg,'Text','u (conditioning value for U1):','Position',[20 100 220 22]);
        edU = uieditfield(dlg,'numeric','Position',[250 100 100 22],'Value',0.5,'Limits',[0.01 0.99]);
        uibutton(dlg,'push','Text','Plot','Position',[140 50 100 26], 'ButtonPushedFcn', @(~,~) conditionalSlice(edU.Value, cols, dlg));
    end

    function conditionalSlice(uval, cols, dlg)
        try
            X = userData.Table{:,cols};
            U = makeUfromMarginals(X, ddlMarg.Value);
            fam = userData.family;
            % pair-specific fit
            switch fam
                case 'Gaussian'
                    R = copulafit('Gaussian', U);
                case 't'
                    [R,nu] = copulafit('t', U);
                otherwise
                    R = []; nu = []; % not needed for Archimedean
            end
            v = linspace(0.001,0.999,400)'; u = max(0.001,min(0.999,uval));
            switch fam
                case 'Gaussian'
                    rho = R(1,2);
                    z1 = norminv(u); mu = rho*z1; s2 = max(1e-8,1-rho^2);
                    z2 = norminv(v);
                    fz = normpdf(z2, mu, sqrt(s2)); fu = normpdf(z2);
                    cond = fz ./ fu;
                case 't'
                    rho = R(1,2); nu0 = max(2,round(nu));
                    z1 = tinv(u, nu0);
                    mu = rho*z1; s2 = (nu0+z1.^2).*(1-rho^2)/(nu0+1);
                    z2 = tinv(v, nu0+1);
                    fz = tpdf((z2-mu)./sqrt(s2), nu0+1) ./ sqrt(s2);
                    fu = tpdf(z2, nu0+1); cond = fz ./ fu;
                otherwise
                    th = copulafit(fam, U);
                    Cline = copulapdf(fam, [repmat(u,numel(v),1) v], th);
                    cond = Cline / trapz(v, Cline);
            end
            figure('Name','Conditional slice U2|U1=u'); plot(v, cond, 'LineWidth',1.4); grid on; xlim([0 1]);
            xlabel(sprintf('v = U_2 | U_1 = %.2f',u)); ylabel('density'); title(sprintf('Conditional density slice (%s)', fam));
            X2 = X(:,2); marg = ddlMarg.Value; pd2 = []; emp=false;
            if ~strcmpi(marg,'Ranks (empirical)')
                try, pd2 = fitdist(X2, lower(marg)); catch, emp=true; end
            else, emp=true; end
            qs = [0.05 0.25 0.5 0.75 0.95];
            if emp, xq = empirical_icdf(X2, qs'); else, xq = icdf(pd2, qs'); end
            hold on; yq = interp1(v, cond, qs, 'linear'); stem(qs, yq, 'filled');
            legend('density','quantile ticks','Location','best');
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Conditional Error');
        end
    end

%% ---------- Simulation ----------
    function simulateFromCopula(src)
        if isempty(userData.U) || (isempty(userData.R) && isempty(userData.theta))
            showAlert(src,'Fit a copula first','Error'); return;
        end
        cols = getSelectedColumns(); d = numel(cols);
        n = max(10, round(edtNSim.Value)); marg = ddlMarg.Value;
        switch lower(userData.family)
            case 't'
                Usim = copularnd('t', userData.R, n, max(2,round(userData.nu)));
            case 'gaussian'
                Usim = copularnd('Gaussian', userData.R, n);
            case 'clayton'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Clayton', userData.theta, n);
            case 'frank'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Frank', userData.theta, n);
            case 'gumbel'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Gumbel', userData.theta, n);
            otherwise
                showAlert(src,'Unknown family','Error'); return;
        end
        X_sim = zeros(n,d);
        for i = 1:d
            x_real = userData.Table{:,cols{i}};
            if strcmpi(marg,'Ranks (empirical)')
                X_sim(:,i) = empirical_icdf(x_real, Usim(:,i));
            else
                pd = fitdist(x_real, lower(marg));
                X_sim(:,i) = icdf(pd, Usim(:,i));
            end
        end
        figure('Name','Simulated Data'); plotmatrix(X_sim); sgtitle(sprintf('Simulated Data from %s copula', userData.family));
    end

    function plotJointPDF(src)
        cols = getSelectedColumns();
        if numel(cols) ~= 2, showAlert(src, 'Select exactly 2 columns for joint PDF plot.', 'Invalid Selection'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','Error'); return; end
        X1 = userData.Table{:,cols{1}}; X2 = userData.Table{:,cols{2}}; marg = ddlMarg.Value;
        U = makeUfromMarginals([X1 X2], marg); U1=U(:,1); U2=U(:,2);
        famKey = userData.family; R2=[]; nu=[]; th=[];
        switch famKey
            case 't'
                [R2,nu]=copulafit('t',[U1 U2]);
            case 'Gaussian'
                R2=copulafit('Gaussian',[U1 U2]);
            case {'Clayton','Frank','Gumbel'}
                th = copulafit(famKey,[U1 U2]);
            otherwise
                showAlert(src,'Unsupported for PDF plot','Error'); return;
        end
        x1 = linspace(min(X1), max(X1), 80); x2 = linspace(min(X2), max(X2), 80);
        [Xgrid, Ygrid] = meshgrid(x1, x2);
        if strcmpi(marg,'Ranks (empirical)')
            U1g = empirical_cdf(X1, Xgrid); U2g = empirical_cdf(X2, Ygrid);
            f1  = empirical_pdf(X1, Xgrid); f2  = empirical_pdf(X2, Ygrid);
        else
            pd1 = fitdist(X1, lower(marg)); pd2 = fitdist(X2, lower(marg));
            U1g = cdf(pd1, Xgrid); U2g = cdf(pd2, Ygrid);
            f1  = pdf(pd1, Xgrid); f2  = pdf(pd2, Ygrid);
        end
        switch famKey
            case 't'
                c = copulapdf('t',[U1g(:) U2g(:)], R2, max(2,round(nu)));
            case 'Gaussian'
                c = copulapdf('Gaussian',[U1g(:) U2g(:)], R2);
            otherwise
                c = copulapdf(famKey,[U1g(:) U2g(:)], th);
        end
        c = reshape(c, size(Xgrid)); jointPDF = c .* f1 .* f2;
        figure('Name','Joint PDF'); surf(x1, x2, jointPDF, 'EdgeColor','none');
        xlabel(cols{1}); ylabel(cols{2}); zlabel('Joint PDF');
        title(sprintf('Joint PDF: %s + %s', userData.family, marg)); view(135,30); colorbar;
    end

   function autoSelectCopula(src)
    cols = getSelectedColumns();
    if isempty(cols) || numel(cols) < 2
        showAlert(src,'Select at least 2 columns','Error'); 
        return; 
    end

    X = userData.Table{:,cols}; 
    marg = ddlMarg.Value; 
    U = makeUfromMarginals(X, marg);

    families = {'Gaussian','t'};
    if size(U,2)==2
        families = [families, {'Clayton','Frank','Gumbel'}];
    end

    results = {};
    for k = 1:numel(families)
        try
            fam = families{k};   % <-- ΔΙΟΡΘΩΣΗ: cell indexing με { }
            if strcmpi(fam,'t')
                [R, nu] = copulafit('t', U);
                logL = sum(log(max(copulapdf('t', U, R, nu), realmin)));
                kparams = numel(R(triu(true(size(R)),1))) + 1;
            elseif strcmpi(fam,'gaussian')
                R = copulafit('Gaussian', U);
                logL = sum(log(max(copulapdf('Gaussian', U, R), realmin)));
                kparams = numel(R(triu(true(size(R)),1)));
            else
                theta = copulafit(fam, U); %#ok<NASGU>
                logL = sum(log(max(copulapdf(fam, U, theta), realmin)));
                kparams = 1;
            end
            nObs = size(U,1);
            AIC = -2*logL + 2*kparams;
            BIC = -2*logL + kparams*log(nObs);
            results(end+1,:) = {fam, logL, AIC, BIC}; %#ok<AGROW>
        catch
            results(end+1,:) = {families{k}, NaN, NaN, NaN}; %#ok<AGROW>
        end
    end

    T = cell2table(results, 'VariableNames', {'Copula','LogL','AIC','BIC'});
    [~, iAIC] = min(T.AIC); 
    [~, iBIC] = min(T.BIC);
    msg = sprintf('Model selection completed:\nBest AIC: %s\nBest BIC: %s', ...
                  T.Copula{iAIC}, T.Copula{iBIC});
    showAlert(src, msg, 'Copula Selection Results');

    figResults = uifigure('Name','Model Selection Table','Position',[320 320 520 220]);
    data = [string(T.Copula), string(round(T.LogL,3)), string(round(T.AIC,3)), string(round(T.BIC,3))];
    uitable(figResults, 'Data', data, 'ColumnName', {'Copula','LogL','AIC','BIC'}, 'Position', [20 20 480 170]);
end

%% ---------- Export / Session ----------
    function exportResults(src)
        if isempty(userData.U), showAlert(src,'Nothing to export (fit first).','Export'); return; end
        [file,path] = uiputfile('copula_results.mat','Save Results'); if isequal(file,0), return; end
        S = struct('U',userData.U,'R',userData.R,'nu',userData.nu,'theta',userData.theta,'family',userData.family,'marg',userData.marg);
        save(fullfile(path,file),'-struct','S'); showAlert(src,'Saved results MAT file.','Export');
    end

    function saveSession(src)
        if isempty(userData.Table), showAlert(src,'No session to save.','Save'); return; end
        [file,path] = uiputfile('copula_session.mat','Save Session'); if isequal(file,0), return; end
        S = struct('Table',userData.Table,'U',userData.U,'R',userData.R,'nu',userData.nu,'theta',userData.theta,'family',userData.family,'marg',userData.marg);
        save(fullfile(path,file),'-struct','S'); showAlert(src,'Session saved.','Save');
    end

    function loadSession(src)
        [file,path] = uigetfile('*.mat','Load Session'); if isequal(file,0), return; end
        S = load(fullfile(path,file));
        if isfield(S,'Table')
            userData.Table = S.Table; names = userData.Table.Properties.VariableNames;
            lstCols.Items = names; if numel(names)>=2, lstCols.Value = names(1:2); else, lstCols.Value=names(1); end
        end
        fields = {'U','R','nu','theta','family','marg'};
        for k=1:numel(fields), if isfield(S,fields{k}), userData.(fields{k}) = S.(fields{k}); end, end
        txtOut.Value = 'Session loaded.';
        if ~isempty(userData.R)
            imagesc(ax,userData.R); colorbar(ax); title(ax,'Copula Correlation (R)');
        end
    end

    function exportAxesPNG(src)
        try
            [file,path] = uiputfile('figure.png','Export current axes'); if isequal(file,0), return; end
            exportgraphics(ax, fullfile(path,file)); showAlert(src,'PNG saved.','Export');
        catch ME
            showAlert(src,ME.message,'Export Error');
        end
    end

    function quickReport(src)
        try
            ts = datestr(now,'yyyymmdd_HHMMSS');
            d = uigetdir(pwd, 'Select folder to create report'); if isequal(d,0), return; end
            rep = fullfile(d, ['report_' ts]); if ~exist(rep,'dir'), mkdir(rep); end
            try, exportgraphics(ax, fullfile(rep,'main_axes.png')); end
            fid = fopen(fullfile(rep,'summary.txt'),'w');
            fprintf(fid,'Copula Toolbox Quick Report\n');
            fprintf(fid,'-----------------------------\n');
            fprintf(fid,'Family: %s\nMarginals: %s\n', userData.family, userData.marg);
            if ~isempty(userData.R)
                fprintf(fid,'R matrix:\n'); fprintf(fid,'%s\n', evalc('disp(userData.R)'));
            end
            if ~isempty(userData.nu), fprintf(fid,'nu: %.4f\n', userData.nu); end
            if ~isempty(userData.theta), fprintf(fid,'theta: %.4f\n', userData.theta); end
            fclose(fid);
            showAlert(src, sprintf('Report created at\n%s', rep), 'Report');
        catch ME
            showAlert(src, ME.message, 'Report Error');
        end
    end

%% ---------- Utilities ----------
    function A = nearestPD(A)
        [V,D] = eig((A+A')/2); D = diag(D); D(D<1e-10) = 1e-10; A = V*diag(D)*V'; A = (A+A')/2;
        s = sqrt(diag(A)); A = A ./ (s*s'); A(1:size(A,1)+1:end) = 1;
    end

    function q = empirical_icdf(x, u)
        x = x(:); u = u(:); [xs,~] = sort(x);
        ranks = (1:numel(x))'/(numel(x)+1);
        q = interp1(ranks, xs, u, 'linear','extrap');
    end

    function Fu = empirical_cdf(x, grid)
        x = x(:); xs = sort(x);
        Fu = arrayfun(@(z) mean(xs<=z, 'omitnan'), grid);
    end

    function fx = empirical_pdf(x, grid)
        x = x(:); if numel(x)<5, fx = zeros(size(grid)); return; end
        bw = 1.06*std(x,'omitnan')*numel(x)^(-1/5);
        if bw<=0 || ~isfinite(bw), bw = max(eps, iqr(x)/1.34* numel(x)^(-1/5)); end
        fx = zeros(size(grid));
        for i=1:numel(x)
            fx = fx + (1/(sqrt(2*pi)*bw)) * exp(-0.5*((grid - x(i))/bw).^2);
        end
        fx = fx / numel(x);
    end

%% ---------- UI helper ----------
    function showAlert(src,msg,title)
        if nargin<3 || isempty(title), title = 'Info'; end
        try
            fig = [];
            if ~isempty(src)
                try, if isvalid(src), fig = ancestor(src,'figure'); end, catch, fig = []; end
            end
            if isempty(fig) || ~isvalid(fig)
                try, if exist('f','var') && isvalid(f), fig = f; end, catch, fig = []; end
            end
            if isempty(fig) || ~isvalid(fig), fig = uifigure('Name',title); end
            uialert(fig, msg, title);
        catch
            try
                warndlg(msg, title);
            catch
                fprintf(2,'%s: %s\n', char(title), char(msg));
            end
        end
    end
end

