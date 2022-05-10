% Damian Strojek 184407
function plot_circ(X, Y, R)
    theta = linspace(0, 2*pi);
    x = R*cos(theta) + X;
    y = R*sin(theta) + Y;
    plot(x,y)
end

% ---------------------------------------------------- 
% Zadanie 1
% ----------------------------------------------------

Edges = sparse([1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7; ...
                4, 6, 4, 5, 3, 6, 5, 7, 6, 5, 6, 4, 4, 7, 6]);

N = 7;              % ilosc stron
d = 0.85;           % parametr przyjety

I = speye(N);
B = sparse(Edges(2,:), Edges(1, :), 1, N, N); % macierz sasiedztwa

L = zeros(1, N);    % same zera
for i = 1:N
    L(i) = sum(B(:,i));
end
A = sparse(diag(1./L));

b = zeros(N, 1);
b(:, 1) = (1 - d)/N;
M = sparse(I - d * B * A);

r = M \ b;  % uklad rownan
bar(r);     % wykres

% ---------------------------------------------------- 
% Zadanie 2
% ----------------------------------------------------

a = 10;             % bok kwadratu
r_max = a/3;        % maksymalny promien
n_max = 100;        % maksymalna ilosc pecherzykow
x = [];y = [];r = [];t = [];

n = 0;                                % aktualna wartosc
rand_count = 0;
rand_counts = [];
area = 0;
areas = [];

while(n <= n_max)
    does_fit = false;
    while(does_fit == false)
        rand_count = rand_count + 1;   
        x_r = rand(1) * a;
        y_r = rand(1) * a;
        r_r = rand(1) * r_max;
         if(x_r + r_r <= a && x_r - r_r > 0 && y_r + r_r <= a && y_r - r_r > 0)
            does_fit = true;
         end
    end

    does_cross = false;
    for i = 1:n
        x_dist = (x_r - x(i));
        y_dist = (y_r - y(i));
        dist = sqrt(x_dist * x_dist + y_dist * y_dist);
        if(dist < r(i) + r_r && dist > abs(r(i) - r_r))
            does_cross = true;
        elseif(dist == r(i) + r_r)
            does_cross = true;
        end
    end

    if(does_cross == false)
        n = n + 1;
        x(n) = x_r;
        y(n) = y_r;
        r(n) = r_r;
        area = area + pi*r_r*r_r;
        areas(n) = area;
         
        % rand_counts(n) = rand_count / n; - takie podejscie nie chcialo
        % mi zadzialac, wiec musialem znalezc alternatywne rozwiazanie
        rand_counts(n) = rand_count;
        rand_count = 0;
        

        fprintf(1, ' %s%5d%s%.3g\r ', 'n =',  n, ' S = ', area);
        pause(0.01);
        axis equal;               % wyrownanie
        plot_circ(x_r, y_r, r_r);
        hold on;
    end
end

figure("Name", "Powierzchnia Całkowita");
semilogx(1:n, areas);
xlabel("Liczba pęcherzyków");
ylabel("Powierzchnia calkowita");
title("Wykres powierzchni calkowitej w zaleznosci od ilosci pecherzykow");
saveas(gcf, "graf1.png");
figure("Name", "Średnia liczba losowań");
loglog(cumsum(rand_counts)./linspace(1, n, n));
xlabel("Liczba pęcherzyków");
ylabel("Średnia liczba losowań");
title("Wykres średniej liczby losowań wielkości pęcherzyków");
saveas(gcf, "graf2.png");
