fx = @(y) cos(y);
fy = @(x) sin(x);
n = 10;
results = fixed_point_contraction(fx, fy, n);
