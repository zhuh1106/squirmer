function test_SLPmatrix
N = 16;
x1c = -5*1i; x2c = 5*1i;
s1.Z = @(t) x1c + 2*cos(t) + 1i*sin(t);
s2.Z = @(t) x2c + 2*cos(t) + 1i*sin(t);
[s1, N, np] = quadr_pan(s1, N, 'p', 'C');
[s2, N, np] = quadr_pan(s2, N, 'p', 'C');

[A1,~] = SLPmatrix(s1,s2,1,0);
A2 = SLPmatrix2(s1,s2);
diff = A1 - A2;
