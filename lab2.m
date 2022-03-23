clc
clear all
close all
diary("log_184407_lab3")

% odpowiednie fragmenty kodu można wykonać poprzez znazaczenie i wciśnięcie F9
% komentowanie/ odkomentowywanie: ctrl+r / ctrl+t


%-----------------
% Zadanie A
%------------------

N = 10;
density = 3; % parametr decydujący o gestosci polaczen miedzy stronami
% tu powinien byc indeks 184407
[Edges] = generate_network(N, density);
%-----------------

% Zadanie B
%------------------
% generacja macierzy I, A, B i wektora b
% macierze A, B i I muszą być przechowywane w formacie sparse (rzadkim)
d = 0.85; 
B = sparse(Edges(2,:), Edges(1,:), 1, N, N);
L = sparse(sum(B));
A = sparse(diag(1./L)); 
I = speye(N);
b = zeros(N, 1);
b(:, 1) = (1-d)/N;

save zadB_184407 A B I b

%-----------------
% Zadanie C
%-----------------
% rozwiazanie ukladu rownan

M = sparse(I - d * B * A);
r = M \ b;

save zadC_184407 r

%-----------------
% Zadanie D
%------------------

clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;
d = 0.85;

for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), 1, N(i), N(i));
    L = sparse(sum(B));
    A = sparse(diag(1./L)); 
    I = speye(N(i));
    b = zeros(N(i), 1);
    b(:, 1) = (1-d)/N(i);
    M = sparse(I - d * B * A);

    tic
    % obliczenia
    r = M \ b;
    czas_Gauss(i) = toc;
end

figure("Name", "Czas Gaussa");
plot(N, czas_Gauss)
title("Zadanie D - Czas rozwiązania metodą Gaussa");
xlabel("Wielkość macierzy N");
ylabel("Czas [s]");
saveas(gcf, "zadD_184407.png");

%------------------
% Zadanie E
%------------------

clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;
d = 0.85;
warunek = 10^(-14);

for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), 1, N(i), N(i));
    L = sparse(sum(B));
    A = sparse(diag(1./L)); 
    I = speye(N(i));
    b = zeros(N(i), 1);
    b(:, 1) = (1-d)/N(i);
    M = sparse(I - d * B * A);
    
    % wzór z pdfa
    r = ones(N(i), 1);
    D = diag(diag(M));
    U = triu(M, 1);
    L = tril(M, -1);
    first = (-D)^(-1) * (L + U);
    second = D^(-1) * b;

    iterations(i) = 1;
    res = 1;
    k = 1;
    tic
    while norm(res) >= warunek
        if N(i) == 1000
            normres(k) = norm(res);
            k = k + 1;
        end
        r = first * r + second;
        res = M * r - b;
        iterations(i) = iterations(i) + 1;
    end
    czas_Jacobi(i) = toc;
end

figure("Name", "Czas Jacobi");
plot(N, czas_Jacobi);
title("Zadanie E - Czas rozwiązania metodą Jacobiego");
xlabel("Wielkość macierzy N");
ylabel("Czas [s]");
saveas(gcf, "zadE_184407_1.png");

figure("Name", "Liczba iteracji Jacobi");
plot(N, iterations);
title("Zadanie E - Liczba iteracji dla metody Jacobiego");
xlabel("Wielkość macierzy N");
ylabel("Liczba iteracji");
saveas(gcf, "zadE_184407_2.png");

figure("Name", "Norma z residuum Jacobi");
semilogy(normres);
title("Zadanie E - Wykres normy residuum Jacobi");
xlabel("Ilość iteracji");
ylabel("Norma residuum")
saveas(gcf, "zadE_184407_3.png");

%------------------
% Zadanie F
%------------------

clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;
d = 0.85;
warunek = 10^(-14);

for i = 1:5
    [Edges] = generate_network(N(i), density);
    B = sparse(Edges(2,:), Edges(1,:), 1, N(i), N(i));
    L = sparse(sum(B));
    A = sparse(diag(1./L)); 
    I = speye(N(i));
    b = zeros(N(i), 1);
    b(:, 1) = (1-d)/N(i);
    M = sparse(I - d * B * A);
    
    % wzór z pdfa
    r = ones(N(i), 1);
    D = diag(diag(M));
    U = triu(M, 1);
    L = tril(M, -1);
    first = -(D + L)^(-1);
    second = (D + L)^(-1) * b;

    iterationsG(i) = 1;
    res = 1;
    k = 1;
    tic
    while norm(res) >= warunek
        if N(i) == 1000
            normres(k) = norm(res);
            k = k + 1;
        end
        r = first * (U*r) + second;
        res = M * r - b;
        iterationsG(i) = iterationsG(i) + 1;
    end
    czas_GaussSeidl(i) = toc;
end

figure("Name", "Czas Gauss-Seidl");
plot(N, czas_GaussSeidl);
title("Zadanie F - Czas rozwiązania metodą Gaussa-Seidla");
xlabel("Wielkość macierzy N");
ylabel("Czas [s]");
saveas(gcf, "zadF_184407_1.png");

figure("Name", "Liczba iteracji Gauss-Seidl");
plot(N, iterationsG);
title("Zadanie F - Liczba iteracji dla metody Gaussa-Seidla");
xlabel("Wielkość macierzy N");
ylabel("Liczba iteracji");
saveas(gcf, "zadF_184407_2.png");

figure("Name", "Norma z residuum Gauss-Siedl");
semilogy(normres);
title("Zadanie F - Wykres normy residuum Gauss-Siedl");
xlabel("Ilość iteracji");
ylabel("Norma residuum");
saveas(gcf, "zadF_184407_3.png");

%------------------
% Zadanie G
%------------------

clc
clear all
close all
load("Dane_Filtr_Dielektryczny_lab3_MN.mat");
r = M \ b;

save zadG_184407_B r

D = diag(diag(M));
U = triu(M, 1);
L = tril(M, -1);
iterationsG1(1) = 0;
iterationsG1(2) = 2;
warunek = 10^(-14);

first = (-D)^(-1) * (L + U);
second = D^(-1) * b;

res = 1;
r = ones(size(D, 1), 1);

while norm(res) >= warunek
    r = first * r + second;
    res = M * r - b;
    iterationsG1(1) = iterationsG1(1) + 1;
    if isnan(norm(res))
        break
    end
end

if isnan(norm(res)) == false
    save zadG_184407_J r
end

%-----------------

first = -(D + L)^(-1);
second = (D + L)^(-1) * b; 

res = 1;
r = ones(size(D, 1), 1);

while norm(res) >= warunek
    r = first * (U*r) + second;
    res = M * r - b;
    iterationsG1(2) = iterationsG1(2) + 1;
    if isnan(norm(res))
        break
    end
end

if isnan(norm(res)) == false
    save zadG_184407_GS r
end


