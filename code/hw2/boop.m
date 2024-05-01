A = readmatrix("output/A.dat");
B = readmatrix("output/B.dat");
A_corner = readmatrix("output/A_corner.dat");
B_corner = readmatrix("output/B_corner.dat");

A - A_corner

B_test = zeros(24);
for i=1:24
    a = 2*(i-1);
    b = 2*i-1;
    a = 2*i - 1;
    b = 2*i;
    if a == 0
        a = 48;
    end
    B_test(:, i) = B_corner(:, a) + B_corner(:, b);
end

B - B_test