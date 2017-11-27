function S = SLPselfmatrix(s) % Stokes SLP via Kress-split Nyst matrix
% Modified based on Laplace case (Barnett 10/18/13). Wu 10/13/14.

% log part
N = numel(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);  % C-# displacements mat, t-s
logpart = -log(abs(d)) + circulant(0.5*log(4*sin([0;s.t(1:end-1)]/2).^2)); % peri log
logpart(diagind(logpart)) = -log(s.sp);                       % diagonal limit
m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2; % Kress Rj(N/2)/4pi
logpart = logpart/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
logpart = logpart .* repmat(s.sp.',[N 1]);  % include speed factors (not 2pi/N weights)
logpart = logpart/2; % Stokes prefactor

% cross product part
r = sqrt(d.*conj(d));    % dist matrix R^{MxN}
d1 = real(d)./r;
d2 = imag(d)./r;
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

t1 = real(s.tang); t2 = imag(s.tang);
cross_part(diagind(cross_part)) = [t1.^2; t2.^2];
cross_part = cross_part(:,[1+end/2:end,1:end/2]);
cross_part(diagind(cross_part)) = [t1.*t2; t2.*t1];
cross_part = cross_part(:,[1+end/2:end,1:end/2]);
    
cross_part = cross_part.*repmat(s.ws(:)', [2*N 2])/4/pi;

S = kron(eye(2),logpart)+cross_part;

function A = circulant(x)
% function A = circulant(x)
%
% return square circulant matrix with first row x
% barnett 2/5/08
x = x(:);
A = toeplitz([x(1); x(end:-1:2)], x);

function i = diagind(A)
% function i = diagind(A)
%
% return diagonal indices of a square matrix, useful for changing a diagonal
% in O(N) effort, rather than O(N^2) if add a matrix to A using matlab diag()
%
% barnett 2/6/08
N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);
