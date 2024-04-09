A = readmatrix("A.dat");
A_old = readmatrix("A_old.dat");

B = readmatrix("B.dat");
B_old = readmatrix("B_old.dat");

x = 1:24;
y = 1:24;

[xm, ym] = meshgrid(x,y);

close all;
figure(1);
contourf(x, y, A - A_old);
colorbar
figure(2);
contourf(x, y, B - B_old);
colorbar
