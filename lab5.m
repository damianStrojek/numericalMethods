function [x, y, f, xp, yp] = lazik(K)
    % K - funkcja zostanie sprobkowana w punktach na plaszczyznie ulozonych w siatke K x K punktow
    % x - wspolrzedne x punktow probkowania
    % y - wspolrzedne y punktow probkowania
    % f - wartosc funkcji w punktach probkowania
    % xp - wspolrzedne x punktow probkowania zgodne z torem jazdy lazika
    % yp - wspolrzedne y punktow probkowania zgodne z torem jazdy lazika

    % generacja polozenia probek
    [x,y] = meshgrid(linspace(0, 1, K),linspace(0, 1, K));

    w1=5; a1=30; x1=0.2; y1=0.3;
    w2=6; a2=40; x2=0.6; y2=1.0;
    w3=7; a3=50; x3=0.9; y3=0.5;
    w4=6; a4=40; x4=0.6; y4=0.1;
    w5=7; a5=70; x5=0.1; y5=0.95;
    w6=-5; a6=10; x6=0.5; y6=0.5;
    f=10+w1*exp(-a1*((x-x1).^2+(y-y1).^2))+w2*exp(-a2*((x-x2).^2+(y-y2).^2))+w3*exp(-a3*((x-x3).^2+(y-y3).^2))+w4*exp(-a4*((x-x4).^2+(y-y4).^2))+w5*exp(-a5*((x-x5).^2+(y-y5).^2))+w6*exp(-a6*((x-x6).^2+(y-y6).^2));

    D = 100; % skala dystansu
    x = x*D;
    y = y*D;

    % utworzenie wektora toru ruchu lazika
    xp = x;
    xp(2:2:end,:) = xp(2:2:end,end:-1:1);
    xp = xp';
    yp = y';

    % macierz -> wektor
    x = reshape(x,K*K,1);
    y = reshape(y,K*K,1);
    f = reshape(f,K*K,1);
    xp = reshape(xp,K*K,1);
    yp = reshape(yp,K*K,1);
end

function [p] = polyfit2d(x, y, f)
    %
    %   A*p = f
    %

    M = sqrt(size(x,1))-1;
    N=M;


    xx = x * ones(1,M+1);
    mm = ones(size(x,1),1)*(0:M);
    yy = y * ones(1,N+1);
    nn = ones(size(y,1),1)*(0:N);

    xm=xx.^mm;
    yn=yy.^nn;

    MN = (M+1)*(N+1);
    A = zeros(MN,MN);

    for k = 1:MN
        
        xmk=xm(k,:);
        ynk=yn(k,:);
        
        A(k,:) = kron(xmk,ynk);
        
    end

    p = A\f;
end

function [FF] = polyval2d(XX, YY, p)
    x=reshape(XX,size(XX,1)*size(XX,2),1);
    y=reshape(YY,size(YY,1)*size(YY,2),1);

    M = sqrt(size(p,1))-1;
    N=M;

    xx = x * ones(1,M+1);
    mm = ones(size(x,1),1)*(0:M);
    yy = y * ones(1,N+1);
    nn = ones(size(y,1),1)*(0:N);

    xm=xx.^mm;
    yn=yy.^nn;

    pr = reshape(p,M+1,N+1);
    f = diag(xm*pr'*yn');

    FF=reshape(f',sqrt(size(f,1)),sqrt(size(f,1)));
end

function [p] = trygfit2d(x, y, f)
    %
    %   A*p = f
    %
    M = sqrt(size(x,1))-1;
    N=M;


    xx = x * ones(1,M+1);
    mm = ones(size(x,1),1)*(0:M);
    yy = y * ones(1,N+1);
    nn = ones(size(y,1),1)*(0:N);

    xm=cos(xx.*mm*pi/max(x));
    yn=cos(yy.*nn*pi/max(y));

    MN = (M+1)*(N+1);
    A = zeros(MN,MN);

    for k = 1:MN
        
        xmk=xm(k,:);
        ynk=yn(k,:);
        
        A(k,:) = kron(xmk,ynk);
        
    end

    p = A\f;
end

function [FF] = trygval2d(XX, YY, p)
    x = reshape(XX,size(XX,1)*size(XX,2),1);
    y = reshape(YY,size(YY,1)*size(YY,2),1);

    M = sqrt(size(p,1))-1;
    N = M;

    xx = x * ones(1,M+1);
    mm = ones(size(x,1),1)*(0:M);
    yy = y * ones(1,N+1);
    nn = ones(size(y,1),1)*(0:N);

    xm = cos(xx.*mm*pi/max(max(XX)));
    yn = cos(yy.*nn*pi/max(max(YY)));

    pr = reshape(p,M+1,N+1);
    f = diag(xm*pr'*yn');

    FF=reshape(f',sqrt(size(f,1)),sqrt(size(f,1)));
end

% Laboratorium 5.1 Damian Strojek 184407

M = 1 + mod(3, 4);
K = [5, 9, 15, 39];

[XX, YY] = meshgrid(linspace(0, 100, 101), linspace(0, 100, 101));
for i = K
    [x, y, f, xp, yp] = lazik(i);

    subplot(2, 2, 1);
    plot(xp, yp, '-o', 'linewidth', 1.25, 'markersize', 2);
    title("Tor ruchu łazika dla K = ", num2str(i));
    ylabel("y [km]");
    xlabel("x [km]");

    subplot(2, 2, 2);
    plot3(x, y, f, 'o', 'markersize', 5);
    title("Wartości próbek dla K = ", num2str(i));
    ylabel("y [km]");
    xlabel("x [km]");
    zlabel("f (x,y)");

    [p] = polyfit2d(x, y, f);
    [FF] = polyval2d(XX, YY, p);
    subplot(2, 2, 3);
    surf(YY, XX, FF);
    shading flat;

    title("Interpolacja wielomianowa dla K = ", num2str(i));
    ylabel("y [km]");
    xlabel("x [km]");
    zlabel("f (x,y)");
    
    [p]  = trygfit2d(x, y, f);
    [FF] = trygval2d(XX, YY, p);
    subplot(2, 2, 4);
    surf(YY, XX, FF);
    shading flat;

    title("Interpolacja trygonometryczna dla K = ", num2str(i));
    ylabel("y [km]");
    xlabel("x [km]");
    zlabel("f (x,y)");

    print(gcf, strcat("K_", num2str(i), ".png"), '-dpng', '-r450');
end

% Laboratorium 5.2 Damian Strojek 184407
               
div_polyval = [];
div_trygval = [];

for i = 5:45
    [XX, YY] = meshgrid(linspace(0, 100, 101), linspace(0, 100, 101));
    [x, y, f, xp, yp] = lazik(i);

    [p] = polyfit2d(x, y, f);
    [FF_p] = polyval2d(XX, YY, p);

    p = trygfit2d(x, y, f);
    [FF_t] = trygval2d(XX, YY, p);

    if i == 5
        FF_p_prev = FF_p;
        FF_t_prev = FF_t;
    else
        div_polyval = [div_polyval, max(max(abs(FF_p - FF_p_prev)))];
        div_trygval = [div_trygval, max(max(abs(FF_t - FF_t_prev)))];
        FF_p_prev = FF_p;
        FF_t_prev = FF_t;
    end
end

plot(6:45, div_polyval);
title("Zbieznosc interpolacji wielomianowej");
ylabel("Maksymalna wartosc roznicy interpolowanych funkcji");
xlabel("Ilosc punktow pomiarowych K");
print (gcf, strcat("Zbieznosc_interpolacja_wielomianowa.png"), '-dpng', '-r450');

plot(6:45, div_trygval);
title("Zbieznosc interpolacji trygonometrycznej");
ylabel("Maksymalna wartosc roznicy interpolowanych funkcji");
xlabel("Ilosc punktow pomiarowych K");
print (gcf, strcat("Zbieznosc_interpolacja_trygonometryczna.png"), '-dpng', '-r450');
