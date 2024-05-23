LHS = complex_readmatrix("output/LHS.mat");
RHS = complex_readmatrix("output/RHS.mat");
A1 = complex_readmatrix("output/A1.mat");
A2 = complex_readmatrix("output/A2.mat");
B1 = complex_readmatrix("output/B1.mat");
B2 = complex_readmatrix("output/B2.mat");

A2aa = A2(1:a, 1:a);
A2ab = A2(1:a, a+1:end);
A2ba = A2(a+1:end, 1:a);
A2bb = A2(a+1:end, a+1:end);
B2aa = B2(1:a, 1:a);
B2ab = B2(1:a, a+1:end);
B2ba = B2(a+1:end, 1:a);
B2bb = B2(a+1:end, a+1:end);

a = 24;
b = 48;
test = zeros(2 * a + b);
test(1:a, 1:a) = A1;
test(1:a, a+1:2*a) = -B1;
test(a+1:end, 1:a) = A2(:, 1:a);
test(a+1:end, a+1:2*a) = B2(:, 1:a);
test(a+1:end, 2*a+1:end) = -B2(:, a+1:end);

sum(LHS - test, 'all')

big_LHS = [[A1, zeros(a), -B1, zeros(a, a+b)];
           [zeros(a) A2aa zeros(a) -B2aa -B2ab];
           [eye(a) -eye(a) zeros(a, 2*a + b)];
           [zeros(a, 2*a) eye(a) eye(a) zeros(a, b)];
           [zeros(b, a), A2ba, zeros(b, a), -B2ba -B2bb]];

big_RHS = [zeros(a, 1); -sum(A2ab, 2); zeros(2*a, 1); -sum(A2bb, 2)];