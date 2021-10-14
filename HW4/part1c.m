clc, clear

N = 100000000;
r = 0:N-1;

tic;
s = 0;
for i = 1:N
    s = s + r(i)*r(i);
end

m = sqrt(s);
time_span = toc;

fprintf("Mag r[%d] is %g\n",N,m);
fprintf("Loop took %g seconds\n", time_span);

tic;
m = sqrt(sum(r.*r));
time_span = toc;
fprintf("Vectors took %g seconds\n", time_span);