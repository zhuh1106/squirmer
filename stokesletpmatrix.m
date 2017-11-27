function A = stokesletpmatrix(t,s) % normal derive of stokeslet matrix
% t = target seg (x,nx cols), s = src seg
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(t.nx, [1 N]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = -real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;