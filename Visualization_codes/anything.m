% MATLAB Code to Save Values of ω and K at Every Iteration for Multiple Flows

clc;
clear;

% Parameters
m = 100;
r = 0;
n = m + r; % Number of flows to draw (adjust this value)
num_iterations = 12; % Number of iterations for each flow

% Initialize random initial conditions for ω and K
%{
omega0_vals = rand(1, n) ; % Random values for ω₀ in [0, 1]
K0_vals = [-power(10,-3*rand(1, m/2+3)+1), power(10,-3*rand(1, m/2-3)+1), -log(2)*ones(1,r)];     % Random values for K₀
%}
%{}
omega0_vals = [rand(1,n-1),0];
K0_vals = [rand(1,n-1), 0.99];
%}
%{
omega0_vals = reshape(ones(n,1)*(0:1/(n-1):1), [1,n*n]).';
K0_vals = reshape(ones(n,1)*(0:1/(n-1):1), [n*n,1]);
%}

% Matrices to store ω and K values
omega_matrix = zeros(n, num_iterations);
K_matrix = zeros(n, num_iterations);

% Create figure
figure;
hold on;
grid on;

% Loop over each flow
K_vec = [];
omega_vec = [];
for j = 1:n
    % Initial values for this flow
    omega_matrix(j, 1) = omega0_vals(j);
    K_matrix(j, 1) = K0_vals(j);

    % Compute recursion
    for i = 2:num_iterations
        %{
        omega_matrix(j, i) = 2 / (omega_matrix(j, i-1) + 1/omega_matrix(j, i-1));
        K_matrix(j, i) = 2 * K_matrix(j, i-1) + 0.5 * ...
            (log(omega_matrix(j, i-1) + 1/omega_matrix(j, i-1)) + log(2));
        %}
        %{}
        omega_matrix(j, i) = sqrt( (K_matrix(j, i-1) + 1/K_matrix(j, i-1) + 2) / (( 1/(omega_matrix(j, i-1)*K_matrix(j, i-1)) + omega_matrix(j, i-1) )*(K_matrix(j, i-1)/omega_matrix(j,i-1) + omega_matrix(j, i-1) )) );
        K_matrix(j, i) = ( power(omega_matrix(j, i-1),2)*K_matrix(j, i-1) + power(K_matrix(j, i-1), 2) ) ...
                    /( power(omega_matrix(j, i-1), 2)*K_matrix(j, i-1) + 1 );
        %}
    end
    K_vec = cat(2,K_vec,K_matrix(j,:));
    omega_vec = cat(2,omega_vec,omega_matrix(j,:));
end

% Plot the flow
%{

[ax,lin2sym_Y,sym2lin_Y] = SLplot(omega_vec,K_vec,'XScale','linear','YScale','symlog');
for j = 1:n
   K = lin2sym_Y(K_matrix(j,:));
    plot(omega_matrix(j, :), K, '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ['Flow ', num2str(j)]);

    % Highlight start point (no legend)
    plot(omega_matrix(j, 1), lin2sym_Y(K_matrix(j, 1)), 'ro', 'MarkerFaceColor', 'r', ...
        'HandleVisibility', 'off');

    % Highlight end point (no legend)
    plot(omega_matrix(j, end), lin2sym_Y(K_matrix(j, end)), 'go', 'MarkerFaceColor', 'g', ...
        'HandleVisibility', 'off');
end

plot([0,1],lin2sym_Y(-log(2)*[1,1]),'--','Color',[.8,.8,.8],'linewidth',2);

ax.XAxis.Limits = [min(omega_vec),max(omega_vec)];

% Customize plot
xlabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 12);
title('Multiple Flows of Recursion on $K_{0}-\omega$ Plane', 'Interpreter', 'latex', 'FontSize', 14);
%legend('show', 'Location', 'best');
%}

%{}
for j = 1:n
   K = K_matrix(j,:);
    plot(omega_matrix(j, :), K, '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', ['Flow ', num2str(j)]);

    % Highlight start point (no legend)
    plot(omega_matrix(j, 1), K_matrix(j, 1), 'ro', 'MarkerFaceColor', 'r', ...
        'HandleVisibility', 'off');

    % Highlight end point (no legend)
    plot(omega_matrix(j, end), K_matrix(j, end), 'go', 'MarkerFaceColor', 'g', ...
        'HandleVisibility', 'off');
end
xlabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 12);
title('Multiple Flows of Recursion on $\alpha-\omega$ Plane', 'Interpreter', 'latex', 'FontSize', 14);
%legend('show', 'Location', 'best')

%}
hold off;
