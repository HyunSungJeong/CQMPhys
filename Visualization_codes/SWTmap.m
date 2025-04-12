[U, J_div_U] = meshgrid(linspace(1, 8, 100), linspace(0, 1/3, 100));
J0 = nan(100,100);
K0 = nan(100,100);
I0 = nan(100,100);
for it1 = 1:100
    for it2 = 1:100
        u = U(it1,it2);
        j = u*J_div_U(it1,it2);
        mu = (u-j)/2;
        if abs(u - 2*j - mu) > 0.01 && abs(u - mu) > 0.01
            j0 = -1/(u-2*j-mu) + 3/(u-mu) + 2/mu;
            k0 = 3/(u-2*j-mu) - 1/(u-mu) + 2/mu;
            i0 = 4*(1/(u-2*j-mu) + 1/(u-mu) + 2/mu);

            if abs(j0) < 15 && abs(k0) < 15
                J0(it1,it2) = j0;
                K0(it1,it2) = k0;
                I0(it1,it2) = i0;
            end

            %{
            if I0(it1,it2) < 0.05
                plot(u,j,'.','Color','blue','MarkerSize',5);
            end
            %}
        end
    end
end

figure;
hold on;
title('$J_{0}$','Interpreter','latex','FontSize',20);
xlim([1,8]);
ylim([0,1/3]);
xlabel('$U$','Interpreter','latex','FontSize',25);
ylabel('$\frac{J}{U}$','Interpreter','latex','FontSize',25);
surf(U, J_div_U, J0, 'EdgeColor', 'none'); % Remove grid lines
view(2); % View from the top (2D)
colormap('jet');
colorbar;
hold off;

figure;
hold on;
title('$K_{0}$','Interpreter','latex','FontSize',20);
xlim([1,8]);
ylim([0,1/3]);
xlabel('$U$','Interpreter','latex','FontSize',25);
ylabel('$\frac{J}{U}$','Interpreter','latex','FontSize',25);
surf(U, J_div_U, K0, 'EdgeColor', 'none'); % Remove grid lines
view(2); % View from the top (2D)
colormap('jet');
colorbar;
hold off;

figure;
hold on;
title('$J_{0} - K_{0}$','Interpreter','latex','FontSize',20);
xlim([1,8]);
ylim([0,1/3]);
xlabel('$U$','Interpreter','latex','FontSize',25);
ylabel('$\frac{J}{U}$','Interpreter','latex','FontSize',25);
surf(U, J_div_U, J0-K0, 'EdgeColor', 'none'); % Remove grid lines
view(2); % View from the top (2D)
colormap('jet');
colorbar;
hold off;

figure;
hold on;
title('$\left|J_{0} - K_{0}\right|$','Interpreter','latex','FontSize',20);
xlim([1,8]);
ylim([0,1/3]);
xlabel('$U$','Interpreter','latex','FontSize',25);
ylabel('$\frac{J}{U}$','Interpreter','latex','FontSize',25);
surf(U, J_div_U, abs(J0-K0), 'EdgeColor', 'none'); % Remove grid lines
view(2); % View from the top (2D)
colormap('jet');
colorbar;
hold off;

figure;
hold on;
title('$\left| I_{0} \right|$','Interpreter','latex','FontSize',20);
xlim([1,8]);
ylim([0,1/3]);
xlabel('$U$','Interpreter','latex','FontSize',25);
ylabel('$\frac{J}{U}$','Interpreter','latex','FontSize',25);
surf(U, J_div_U, abs(I0), 'EdgeColor', 'none'); % Remove grid lines
view(2); % View from the top (2D)
colormap('jet');
colorbar;
hold off;