clear;

N = 20;
Lambda = 4;
X = power(sqrt(Lambda), 0:-0.01:-N );
Y = zeros(1,numel(X));
alpha = 0.7;

for it = 1:N
    omega_0 = 100*power(alpha*sqrt(Lambda), -it);
    for itX = 1:numel(X) 
        Y(itX) = Y(itX) + 1/(omega_0 * sqrt((omega_0^2 - X(itX)^2)^2 + omega_0^2 * 0.000001));
    end
end

figure;
hold on;
set(gca,'XScale','log','YScale','log');
plot(X,Y);
hold off;