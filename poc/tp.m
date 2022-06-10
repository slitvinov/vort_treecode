x = 0:0.52360:5.7596;
y = [408 89 -66 10 338 807 1238 1511 1583 1462 1183 804];

plot(x, y, 'o');

m = length(y);
n = floor((m+1)/2);

z = fft(y)/m;

a0 = z(1);
an = 2*real(z(2:n));
bn = -2*imag(z(2:n));
al = z(n+1);

px = 0:0.01:2*pi;
k = 1:length(an);
py = a0 + an*cos(k'*px) ...
	+ bn*sin(k'*px) ...
	+ al*cos(n*px);

plot(x, y, 'o', px, py, '-')
