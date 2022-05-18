function [wsp_wielomianu, x_approx] = aproksymacjaWiel(n, x, N)

  mu = [mean(n) std(n)];
  n2 = (n-mu(1))/mu(2);
  wsp_wielomianu = polyfit(n2,x,N);
  x_approx = polyval(wsp_wielomianu,n2);

end

function [x_aprox] = aprox_tryg(n,x,N)

  n = n*(pi/max(n));
  S = zeros(N+1,N+1);

  % generacja macierzy S
  for i = 1:N+1
      for j = 1:N+1
          S(i,j) = sum(cos((i-1).*n).*cos((j-1).*n));
      end
  end

  t = zeros(N+1, 1);
  for k = 1:N+1
      t(k,1) = sum(x.*cos((k-1)*n));
  end


  % Rozwiazanie ukladu rownan Sc = t
  c = S \ t;

  c1 = cos((0:N)'*n);

  x_aprox = (c1' * c).';

end

% ---------------------------------------
% maindron.m
% ---------------------------------------
clc
clear all
close all

warning('off','all')

load trajektoria1

N = 90;
[ wsp_wielomianu, xa ] = aproksymacjaWiel(n,x,N);  % aproksymacja wsp. 'x'.

% ---------------------------------------
% zadanie2.m
% ---------------------------------------   
clc
clear all
close all

M = 1 + mod(3,4);

load trajektoria1

plot3(x,y,z,'o');
grid on
axis equal
title("Trajektoria rzeczywista");
xlabel("x axis [m]");
ylabel("y axis [m]");
zlabel("z axis [m]");

print(gcf, '184407_strojek_zad2.png', '-dpng', '-r450')
             
% ---------------------------------------
% zadanie3.txt
% ---------------------------------------
% Aproksymacja jest w tym przypadku lepsza ponieważ dane dotyczące położenia drona są obarczone błędem. 
% Gdybyśmy użyli metody interpolacji wynik byłby dopasowany do węzłów,
% które grubo odbiegają od rzeczywistych wartości przez co wynik byłby przekłamany.
             
% ---------------------------------------
% zadanie4.m
% ---------------------------------------
clc
clear all
close all

load trajektoria1

plot3(x, y, z, 'o');
grid on
axis equal
hold on

N = 50;
[~, ax] = aproksymacjaWiel(n, x, N);
[~, ay] = aproksymacjaWiel(n, y, N);
[~, az] = aproksymacjaWiel(n, z, N);

plot3(ax, ay, az, 'lineWidth', 4);
title("Trajektoria rzeczywista i aproksymowana dla N=50");
xlabel("x axis [m]");
ylabel("y axis [m]");
zlabel("z axis [m]");
hold off
print (gcf, '184407_strojek_zad4.png', '-dpng', '-r450')
             
% ---------------------------------------
% zadanie5.m
% ---------------------------------------
clc
clear all
close all

load trajektoria2

plot3(x, y, z, 'o');
grid on
axis equal
hold on

N = 60;
[~, ax] = aproksymacjaWiel(n, x, N);
[~, ay] = aproksymacjaWiel(n, y, N);
[~, az] = aproksymacjaWiel(n, z, N);

plot3(ax, ay, az, 'lineWidth', 4);

title("Trajektoria rzeczywista i aproksymowana dla N=60");
xlabel("x axis [m]");
ylabel("y axis [m]");
zlabel("z axis [m]");
hold off

print (gcf, '184407_strojek_zad5.png', '-dpng', '-r450')

err = [];
M = size(n, 2);
for N = 1:71
    [~, ax] = aproksymacjaWiel(n, x, N);
    [~, ay] = aproksymacjaWiel(n, y, N);
    [~, az] = aproksymacjaWiel(n, z, N);
    
    err = [err, (sqrt(sum((x - ax).^2))) / M + (sqrt(sum((y - ay).^2)) / M) + (sqrt(sum((z - az).^2)) / M)];
end

figure;
semilogy(err);
grid on;
title("Błąd aproksymacji wielomianowej");
xlabel("N");
ylabel("Wartość błędu");
print (gcf, '184407_strojek_zad5_b.png', '-dpng', '-r450')

% ---------------------------------------
% zadanie7.m
% ---------------------------------------
clc
clear all
close all

load trajektoria2

figure;
subplot(1, 2, 1);
plot3(x, y, z, 'o');
grid on
axis equal
hold on

N = 60;
ax = aprox_tryg(n, x, N);
ay = aprox_tryg(n, y, N);
az = aprox_tryg(n, z, N);

plot3(ax, ay, az, 'lineWidth', 4);
title("N=60");
xlabel("x axis [m]");
ylabel("y axis [m]");
zlabel("z axis [m]");
hold off

subplot(1, 2, 2)
plot3(x, y, z, 'o')
grid on
axis equal
hold on

M = size(n, 2);
epsilon = 10^-3;
i = 0;

while true
  i = i + 1;
  ax = aprox_tryg(n, x, i);
  ay = aprox_tryg(n, y, i);
  az = aprox_tryg(n, z, i);
    
  err = (sqrt(sum((x - ax).^2)) / M) 
        + (sqrt(sum((y - ay).^2)) / M) 
        + (sqrt(sum((z - az).^2)) / M);
  
  if err <= epsilon
    break
  end
end

ax = aprox_tryg(n, x, i);
ay = aprox_tryg(n, y, i);
az = aprox_tryg(n, z, i);

plot3(ax, ay, az, 'lineWidth', 4);
title("N=auto");
xlabel("x axis [m]");
ylabel("y axis [m]");
zlabel("z axis [m]");
hold off

print (gcf, '184407_strojek_zad7.png', '-dpng', '-r450')

err = [];
for N = 1:71
  ax = aprox_tryg(n, x, N);
  ay = aprox_tryg(n, y, N);
  az = aprox_tryg(n, z, N);
    
  err = [err, (sqrt(sum((x - ax).^2)) / M) + (sqrt(sum((y - ay).^2)) / M) + (sqrt(sum((z - az).^2)) / M)];
end

figure;
semilogy(err);
grid on;
title("Błąd aproksymacji trygonometrycznej");
xlabel("N");
ylabel("Wartość błędu");
print (gcf, '184407_strojek_zad7_b.png', '-dpng', '-r450')
