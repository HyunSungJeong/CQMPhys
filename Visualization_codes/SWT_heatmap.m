% Parameters
U = 10;
Hyb = 1;

% Define grid (extended to include negative J, JL)
J_vals  = linspace(-U, U, 1000);      % explore negative and positive J
JL_vals = linspace(-U, U, 1000);      % explore negative and positive J_L
[J, JL] = meshgrid(J_vals, JL_vals);

% Compute candidate bounds for mu
expr1 = U - 2*J + JL/4;
expr2 = U - (3/4)*JL;
expr3 = U + JL/4;

minExpr = min(min(expr1, expr2), expr3);

% Define mu from condition in the problem
mu = 0.5 * (minExpr - JL/4);
gap = 0.5 * (minExpr + JL/4);

% Check inequality alpha:  -JL/4 < mu < minExpr
cond1 = (-JL/4 < mu);
cond2 = (mu < minExpr);
validMask = cond1 & cond2;

% Define I
I = 1 ./ (U - 2*J + JL/4 - mu) + ...
    1 ./ (U - 3*JL/4 - mu) + ...
    2 ./ (JL/4 + mu);

J0 = -1 ./ (U-2*J + JL/4 - mu) ...
    +2 ./ (U - 3*JL/4 - mu) ...
    +1 ./ (U + JL/4 - mu) ...
    +2 ./ (JL/4 + mu);
J0 = J0 ./ 2;

K_perp = 3 ./ (U-2*J + JL/4 - mu) ...
         -1 ./ (U + JL/4 - mu) ...
         +2 ./ (JL/4 + mu);
K_perp = K_perp ./ 2;

K_z = 3 ./ (U-2*J + JL/4 - mu) ...
      -2 ./ (U - 3*JL/4 - mu) ...
      +1 ./ (U + JL/4 - mu) ...
      +2 ./ (JL/4 + mu);
K_z = K_z ./ 2;

I = Hyb*I;  J0 = Hyb*J0;    K_perp = Hyb*K_perp;   K_z = Hyb*K_z;

plotLOGI = true;
plotJ = false;
plotKperp = false;
plotKz = false;
plotmu = false;
plotgap = false;
plot_JdivK = true;
plotLOGratio = false;
plotratio = true;

%% Heatmap of log10(I0/K_perp)
if plotLOGI
    I_plot = log10(I ./ K_perp);
    I_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    minVal = min(I_plot(validMask), [], 'all');
    
    % Clip values above 10
    I_plot(I_plot > 1) = 1;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, I_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; jet(256)]); % black for invalid, parula otherwise
    clim([minVal-0.1, 0.7]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$ \log_{10}(I_{0} / K_{\perp}) \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% Heatmap of J0

if plotJ
    J_plot = J0;
    J_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    minVal = min(J_plot(validMask), [], 'all');
    
    % Clip values above 10
    J_plot(J_plot > 0.6) = 0.6;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, J_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]); % black for invalid, parula otherwise
    clim([minVal-0.3, 1.1]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$J_{0} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% Heatmap of K_perp

if plotKperp
    Kperp_plot = K_perp;
    Kperp_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    %
    minVal = min(Kperp_plot(validMask), [], 'all');
    
    % Clip values above 10
    Kperp_plot(Kperp_plot > 0.6) = 0.6;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, Kperp_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]); % black for invalid, parula otherwise
    clim([minVal-0.3, 1.1]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$K_{\perp} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% Heatmap of K_z

if plotKz
    Kz_plot = K_z;
    Kz_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    %
    minVal = min(Kz_plot(validMask), [], 'all');
    
    % Clip values above 10
    Kz_plot(Kz_plot > 0.6) = 0.6;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, Kz_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]); % black for invalid, parula otherwise
    clim([minVal-0.3, 1.1]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$K_{z} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% Heatmap of J0/K_perp, J0/K_z

if plot_JdivK

    % plot J0/K_perp

    JdivK_plot = J0 ./ K_perp;
    JdivK_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    minVal = min(JdivK_plot(validMask), [], 'all');
    maxVal = max(JdivK_plot(validMask), [], 'all');
    
    % Clip values above 10
    JdivK_plot(JdivK_plot > 0.6) = 0.6;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, JdivK_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; jet(256)]); % black for invalid, parula otherwise
    clim([minVal, 0.7]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$J_{0} / K_{\perp} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);


    % plot J0/K_z

    JdivK_plot = J0 ./ K_z;
    JdivK_plot(~validMask) = NaN;
    
    % Determine min of log10|I| over valid region
    minVal = min(JdivK_plot(validMask), [], 'all');
    maxVal = max(JdivK_plot(validMask), [], 'all');
    
    % Clip values above 10
    JdivK_plot(JdivK_plot > 0.6) = 0.6;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, JdivK_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; jet(256)]); % black for invalid, parula otherwise
    clim([minVal, 0.7]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$J_{0} / K_{z} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% --- Heatmap of mu ---

if plotmu
    mu_plot = mu;
    mu_plot(~validMask) = NaN;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, mu_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]);
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$\mu \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% --- Heatmap of gap/sqrt(Hyb) ---

if plotgap
    gap_plot = gap/sqrt(Hyb);
    gap_plot(~validMask) = NaN;
    
    minVal = min(gap(validMask), [], 'all');
    
    figure;
    imagesc(J_vals/U, JL_vals/U, gap_plot);
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]);
    %clim([minVal, 2.5]);
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$\mathrm{gap} \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% --- Find zero/minimum of ratio in valid region ---
I_valid = I;
I_valid(~validMask) = NaN;
J0_valid = J0;
J0_valid(~validMask) = NaN;

ratio = nan(size(I_valid));

for it1 = 1:size(I_valid, 1)
    for it2 = 1:size(I_valid, 2)
        MaxCoupling = max(abs([J0(it1,it2), K_perp(it1,it2), K_z(it1,it2)]));

        if isinf(MaxCoupling)
            MaxCoupling = 1e2;
            ratio(it1,it2) = NaN;
        else
            ratio(it1,it2) = I_valid(it1,it2) / MaxCoupling;
        end

    end
end

%ratio = abs(I_valid(:) ./ K_perp(:) );
%[~, idx] = min(abs(I_valid(:)));
[mindfd, idx] = min(ratio(:));

J_opt  = J(idx);
JL_opt = JL(idx);

disp(mu(idx)+JL(idx)/4);
disp(minExpr(idx) - mu(idx));

fprintf('Approximate zero/minimum of I at: J = %.4f, J_L = %.4f\n', J_opt, JL_opt);
fprintf('I(J,J_L) = %.4e, J0(J,J_L) = %.4e, K_p(J,J_L) = %.4e, K_z(J,J_L) = %.4e, mu = %.4f\n', I(idx), J0(idx), K_perp(idx), K_z(idx), mu(idx));

%% --- Heatmap of log10(ratio) ---

if plotLOGratio
    % Determine min of log10|I| over valid region
    %minVal = min(log10(ratio(validMask)), [], 'all');
    minVal = -1.5;
    
    figure;
    imagesc(J_vals/U, JL_vals/U, log10(ratio));
    axis xy;
    colorbar;
    colormap([0 0 0; parula(256)]); % black for invalid, parula otherwise
    clim([minVal, 1]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$\log_{10}(\mathrm{ratio}) \ (\mathrm{black = invalid \ region})$', 'Interpreter', 'latex', 'FontSize', 18);
end

%% --- Heatmap of ratio ---

if plotratio
    % Determine min of log10|I| over valid region
    %minVal = min(ratio(validMask), [], 'all');
    minVal = min(ratio(validMask), [], 'all');
    maxVal = max(ratio(validMask), [], 'all');
    
    figure;
    imagesc(J_vals/U, JL_vals/U, ratio);
    axis xy;
    colorbar;
    colormap([0 0 0; jet(256)]); % black for invalid, parula otherwise
    clim([minVal, maxVal]);              % scale between min and 1
    xlabel('$J/U$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$J_{\mathrm{L}}/U$', 'Interpreter', 'latex', 'FontSize', 18);
    title('$\mathrm{mean}(I_{\perp}, I_{z}) / \max(J_{0}, K_{\perp}, K_{z})$', 'Interpreter', 'latex', 'FontSize', 18);
end