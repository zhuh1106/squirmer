function [A,T] = SLPmatrix(t,s,mu,a) % single-layer 2D Stokes kernel vel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow component
% t = target seg (x cols), s = src seg, a = optional translation of src seg
% No option for self-int, gives Inf on diag.   3/2/14
% 2nd output is traction matrix, needs t.nx normal (C-#); no self-int either.

if nargin==0&&nargout==0
    STR=dbstack;
    STR=input(['Do you want to test ',STR.name,'.m? Use y/n: '],'s');
    if STR=='y'
        testcode
        return
    end
end

if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); d2 = imag(r);
Ilogr = -log(abs(r));  % log(1/r) diag block
c = 1/(4*pi*mu);       % factor from Hsiao-Wendland book, Ladyzhenskaya
A12 = c*d1.*d2.*irr;   % off diag vel block
A = [c*(Ilogr + d1.^2.*irr), A12;                         % u_x
     A12,                c*(Ilogr + d2.^2.*irr)];         % u_y
A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);              % quadr wei
if nargout>1           % traction (negative of DLP vel matrix w/ nx,ny swapped)
  rdotn = d1.*repmat(real(t.nx), [1 N]) + d2.*repmat(imag(t.nx), [1 N]);
  rdotnir4 = rdotn.*(irr.*irr); clear rdotn
  A12 = -(1/pi)*d1.*d2.*rdotnir4;
  T = [-(1/pi)*d1.^2.*rdotnir4,   A12;                   % own derivation
     A12,                      -(1/pi)*d2.^2.*rdotnir4];
  T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end
end

function testcode
mu=0.7; %viscosity
a=linspace(0,2*pi,257)';
a=a(2:end);
s.x=[cos(a) sin(a)+cos(2*a)/2]*[1;1i]; % vesicle position
s=quadr(s); % Computes properties of the vesicle
t.x=rand(10,2)*[1;1i]+3*1i; % target points
tau=[exp(sin(a)),exp(cos(a))]*[1;1i]; % density
TAU=[real(tau);imag(tau)];% density as a real vector

% Test 1: Comparing far interactions with StokesDcloseeval
z1=SLPmatrix(t,s,mu)*TAU;
z1=z1(1:end/2)+1i*z1(end/2+1:end);
z2=StokesScloseeval(t.x,s,tau,'e')/mu; % close evaluation scheme
error1=max(abs(z1-z2));
if error1<1e-14
    disp('Test 1: Success')
else
    disp('Test 1: Failure')
end
end