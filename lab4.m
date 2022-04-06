% Damian Strojek 184407
clc
close all

% ------------------------------------------------------------------------

clear
a = 0;   
b = 50;

[xvect, xdif, ~, it_cnt] = bisect(a, b, 1e-12, @compute_impedance);

semilogy(1:it_cnt, xvect);
title("1.1");
ylabel("Przybliżona wartość prędkości kątowej omega [rad/s]");
xlabel("Numer iteracji i");
saveas(gcf, 'zad1_1_184407_przyblizenia_bisekcja.png');

plot(1:it_cnt, xdif);
title("1.2");
ylabel("Różnca pomiędzy obecną i poprzednią wartością prędkości kątowej");
xlabel("Numer iteracji i");
saveas(gcf, "zad1_2_184407_roznice_bisekcja.png");

[xvect, xdif, ~, it_cnt] = secant(a, b, 1e-12, @compute_impedance);

semilogy(1:it_cnt, xvect);
title("1.3");
ylabel("Przybliżona wartość prędkości kątowej omega [rad/s]");
xlabel("Numer iteracji i");
saveas(gcf, 'zad1_3_184407_przyblizenia_sieczne.png');

plot(1:it_cnt, xdif);
title("1.4");
ylabel("Różnca pomiędzy obecną i poprzednią wartością prędkości kątowej");
xlabel("Numer iteracji i");
saveas(gcf, "zad1_4_184407_roznice_sieczne.png");

% ------------------------------------------------------------------------

clear
a = 1;
b = 60000;

[xvect, xdif, ~, it_cnt] = bisect(a, b, 1e-3, @time);

semilogy(1:it_cnt, xvect);
title("2.1 - wartość kolejnego przybliżenia x_i w zależności od numeru iteracji i przy użyciu metody bisekcji");
ylabel("Liczba parametrów wejściowych");
xlabel("Numer iteracji i");
saveas(gcf, "zad2_1_184407_przyblizenia_bisekcja.png");

plot(1:it_cnt, xdif);
title("2.2 - różnice pomiędzy wartościami x w kolejnych iteracjach przy użyciu metody bisekcji");
ylabel("Różnica pomiędzy obecną i poprzednią wartością");
xlabel("Numer iteracji i");
saveas(gcf, "zad2_2_184407_roznice_bisekcja.png");

[xvect, xdif, ~, it_cnt] = secant(a, b, 1e-3, @time);

semilogy(1:it_cnt, xvect);
title("2.3 - wartość kolejnego przybliżenia x_i w zależności od numeru iteracji i przy użyciu metody siecznych");
ylabel("Liczba parametrów wejściowych");
xlabel("Numer iteracji i");
saveas(gcf, "zad2_3_184407_przyblizenia_sieczne.png");

plot(1:it_cnt, xdif);
title("2.4 - różnice pomiędzy wartościami x w kolejnych iteracjach przy użyciu metody siecznych");
ylabel("Różnica pomiędzy obecną i poprzednią wartością");
xlabel("Numer iteracji i");
saveas(gcf, "zad2_4_184407_roznice_sieczne.png");

% ------------------------------------------------------------------------

clear
a = 0;
b = 50;

[xvect, xdif, ~, it_cnt] = bisect(a, b, 1e-12, @speed);

semilogy(1:it_cnt, xvect);
title("3.1 - wartość kolejnego przybliżenia t w zależności od numeru iteracji przy użyciu metody bisekcji");
ylabel("Przybliżona wartość czasu [s]");
xlabel("Numer iteracji i");
saveas(gcf, "zad3_1_184407_przyblizenia_bisekcja.png");

plot(1:it_cnt, xdif);
title("3.2 - różnice pomiędzy wartościami t w kolejnych iteracjach przy użyciu metody bisekcji");
ylabel("Przybliżona wartość czasu [s]");
xlabel("Numer iteracji i");
saveas(gcf, "zad3_2_184407_roznice_bisekcja.png");

[xvect, xdif, ~, it_cnt] = secant(a, b, 1e-12, @speed);

semilogy(1:it_cnt, xvect);
title("3.3 - wartość kolejnego przybliżenia t w zależności od numeru iteracji przy użyciu metody siecznych");
ylabel("Przybliżona wartość czasu [s]");
xlabel("Numer iteracji i");
saveas(gcf, "zad3_3_184407_przyblizenia_sieczne.png");

plot(1:it_cnt, xdif);
title("3.4 - różnice pomiędzy wartościami t w kolejnych iteracjach przy użyciu metody siecznych");
ylabel("Przybliżona wartość czasu [s]");
xlabel("Numer iteracji i");
saveas(gcf, "zad3_4_184407_roznice_sieczne.png");

options = optimset('Display', 'iter');
fzero(@tan, 6, options);
fzero(@tan, 4.5, options);

function difference = compute_impedance(omega)
    R = 725;
    C = 8*1e-5;
    L = 2;
    Z = 1/(sqrt((1/R^2) + (omega*C - 1/(omega*L))^2));
    difference = Z - 75;
end

function difference = speed(t)
    g = 9.81;
    q = 2700;
    m = 150000;
    u = 2000;
    v = u * log(m/(m-q*t)) - g*t;
    difference = v - 750;
end

function difference = time(N)
    t = (N^1.43 + N^1.14)/1000;
    difference = t - 5000;
end

function [xvect, xdif, f_x, it_cnt] = bisect(a, b, eps, fun)
    x_temp = a;
    xvect = [];
    xdif = [];
    f_x = [];

    for i = 1:1000
        x = (a + b)/2;  
        % bisection algorithm 
        f_x_temp = feval(fun, x);
        xvect(i) = x;
        xdif(i) = abs(x_temp - x);
        f_x(i) = f_x_temp;
        x_temp = x;
        if abs(f_x_temp) <= eps || abs(a - b) > eps
            it_cnt = i;
            return;
        elseif f_x_temp * feval(fun, a) < 0
            b = x;
        else   
            a = x;
        end
    end
end

function [xvect, xdif, f_x, it_cnt] = secant(a, b, eps, fun)
    x_1 = a;
    x_2 = b;

    xvect = [];
    xdif = [];
    f_x = [];

    for i = 1:1000
        f_x_temp = feval(fun, x_2);
        % derivation of the method
        x_3 = x_2 - (f_x_temp * (x_2-x_1) / (f_x_temp - feval(fun, x_1)));
        xvect(i) = x_3;
        xdif(i) = abs(x_3 - x_2);
        f_x(i) = feval(fun, x_3);

        if abs(f_x(i)) < eps
            it_cnt = i;
            return;
        end

        x_1 = x_2;
        x_2 = x_3;
    end
end
