function squirmer_periodic_main2
% for the purpose of verifying velocity add up
clear all
close all

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

N = 400; %80; % pts per top and bottom wall (enough for 1e-15 in this case)
U.Z = @(t) t + 1i*(1.5+sin(t)); U.Zp = @(t) 1 + 1i*cos(t);
U.Zpp = @(t) -1i*sin(t); U = quadr(U,N); % t=0 is left end, t=2pi right end
D.Z = @(t) t + 1i*(-1.5+cos(2*t)); D.Zp = @(t) 1 - 2i*sin(2*t);
D.Zpp = @(t) - 4i*cos(2*t); D = quadr(D,N); % same direction (left to right)
U.nx = -U.nx; U.cur = -U.cur; % correct for sense of U, opp from periodicdirpipe

for i = 3:3

N1 = i*100;
N2 = i*100;
N3 = i*100;
dt = 0.5;
theta1 = 0; theta2 = 0;
beta = -0.4;
x1c = 2;

% to hold all squirmer
s = {};
% define single particle geometry attributes
s{1}.Z  = @(t) x1c + 1/2*cos(t) + 1/4*1i*sin(t);
s{1}.Zr = @(t) real(s{1}.Z(t) - x1c);
s{1}.Zi = @(t) imag(s{1}.Z(t) - x1c);
s{1}.xc = x1c;
[s{1}, N1, np] = quadr_pan(s{1}, N1, 'p', 'C');
s{1}.n = N1;
s{1}.theta = 0;



numOfP = 1; %num of squirmer
numOfQ = N1; %num of quadrature
rhs = zeros(2*numOfQ+3*numOfP,1);

%     plot(real(s{1}.x),imag(s{1}.x))
%     hold on
%     axis equal
%     xlim([-5 20])
%     ylim([-8 8])
numOfP = numel(s);
uInf = {};
for i = 1:numOfP
    uInf{i}.ug = zeros(2*s{i}.n,1); 
end
[s,mu] = squirmer_solver(s, beta ,dt, uInf);

% figure(),
% plot(U.x,'k.-'); hold on;  plot(D.x,'k.-'); 
% plot(real(s{1}.x),imag(s{1}.x),'b'); 

nx = 60; gx = 2*pi*((1:nx)-0.5)/nx; ny = nx; gy = gx - pi; % plotting grid
[xx yy] = meshgrid(gx,gy); t.x = xx(:)+1i*yy(:); Mt = numel(t.x);
ug = SLPmatrix(t,s{1},1) * mu; 
ug = ug + SLPmatrix(t,s{1},1,-2*pi) * mu; 
ug = ug + SLPmatrix(t,s{1},1,2*pi) * mu; 
u1 = reshape(ug(1:Mt),size(xx)); u2 = reshape(ug(Mt+(1:Mt)),size(xx));
    
[s_test,mu] = squirmer_solver(s, beta ,dt, uInf);   


Smu = {};
Smu{1}.mu = mu;
[u1c,u2c,di]=stokesperivelpipe(Smu, s);


xForPlot = s_test{1}.x;
[IN, ON] = inpolygon(real(t.x),imag(t.x),real(xForPlot),imag(xForPlot));
ii = ~IN;
ii = reshape(ii,size(xx));

figure()
plot(real(s_test{1}.x),imag(s_test{1}.x),'g'); hold on;
plot(U.x,'k.-'); hold on;  plot(D.x,'k.-'); 
quiver(gx,gy, (u1c+u1).*di.*ii,(u2c+u2).*di.*ii, 3);
% uInf = stokesperivelpipe2(Smu, s);
% [s,mu] = squirmer_solver(s, beta ,dt, uInf);
% plot(real(s{1}.x),imag(s{1}.x),'r'); hold off;
axis([0 2*pi -pi pi]);
% pos(:,i) = [s{1}.xc,s{2}.xc,s{3}.xc];
end
