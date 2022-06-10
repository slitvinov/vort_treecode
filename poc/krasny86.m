1;
function [u, v] = rhs(dx, dy)
  global delta
  D = cosh(2 * pi * dy) .- cos(2 * pi * dx) + delta^2;
  u = -sinh(2 * pi * dy) ./ D;
  v = sin(2 * pi * dx) ./ D;
endfunction

function dr = fun(t, r)
  global N
  x = r(1:N);
  y = r(N+1:end);
  dr = zeros(2 * N, 1);
  for j=1:N
    [u v] = rhs(x(j) - x(1:end ~= j),
		y(j) - y(1:end ~= j));
    dr(j) = sum(u);
    dr(j + N) = sum(v);
  endfor
  dr /= 2 * N;
endfunction

global N delta
delta = 0.5;
N = 400;
dG = 1 / N;
j = (1:N)';
G = dG * (j - 1);
x = G + 0.01 * sin(2 * pi * G);
y = -0.01 * sin(2 * pi * G);
[t R] = ode45(@fun, [0 4], [x; y]);

plot(R(1, 1:N), R(1, N+1:end), 'b-', R(end, 1:N), R(end, N+1:end), 'ro-');
axis([0 1 -0.275 0.275])
pause()
