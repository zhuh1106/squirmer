function three_particle2

clear all
close all

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

qtype = 'p';
kappa = 0;
amp = 20;
numOfP = 1; %num of squirmer

size = 4;
N1 = 64*size;
N2 = 64*size;
N3 = 64*size;
dt = 0.05;
theta1 = 0; theta2 = 0;
x1c = -5*1i; x2c = 5*1i; x3c = 0;

% to hold all squirmer
s = {};
% define single particle geometry attributes
s{1}.Z  = @(t) x1c + 2*cos(t) + 1i*sin(t);
s{1}.Zr = @(t) real(s{1}.Z(t) - x1c);
s{1}.Zi = @(t) imag(s{1}.Z(t) - x1c);
s{1}.xc = x1c;
% [s{1}, N1, np] = quadr_pan(s{1}, N1, 'p', 'C');
if qtype=='g'
    s{1} = quadr(s{1}, N1);
else
    [s{1}, N1, np] = quadr_pan(s{1}, N1, 'p', 'C');
end
s{1}.n = N1;
s{1}.theta = 0;

s{2}.Z  = @(t) x2c + 2*cos(t) + 1i*sin(t);
s{2}.Zr = @(t) real(s{2}.Z(t) - x2c);
s{2}.Zi = @(t) imag(s{2}.Z(t) - x2c);
s{2}.xc = x2c;
% [s{2}, N2, np] = quadr_pan(s{2}, N2, 'p', 'C');
if qtype=='g'
    s{2} = quadr(s{2}, N2);
else
    [s{2}, N2, np] = quadr_pan(s{2}, N2, 'p', 'C');
end
s{2}.n = N2;
s{2}.theta = 0;

s{3}.Z  = @(t) x3c + 2*cos(t) + 1i*sin(t);
s{3}.Zr = @(t) real(s{3}.Z(t) - x3c);
s{3}.Zi = @(t) imag(s{3}.Z(t) - x3c);
s{3}.xc = x3c;
% [s{3}, N3, np] = quadr_pan(s{3}, N3, 'p', 'C');
if qtype=='g'
    s{3} = quadr(s{3}, N3);
else
    [s{3}, N3, np] = quadr_pan(s{3}, N3, 'p', 'C');
end
s{3}.n = N3;
s{3}.theta = 0;

uInf = {};
uInf{1}.ug = 0; uInf{2}.ug = 0; uInf{3}.ug = 0;

% rhs condition

numOfQ = N1; %num of quadrature
% rhs = zeros(2*numOfQ+3*numOfP,1);
beta = -0.4;
% beta = 0.3;

% make movie
v = VideoWriter('swimmer.avi');
open(v);
    plot(real(s{1}.x),imag(s{1}.x))
    hold on, plot(real(s{2}.x),imag(s{2}.x))
    plot(real(s{3}.x),imag(s{3}.x))
    hold off
    axis equal
    xlim([-5 20])
    ylim([-8 8])
    pause(0.05)
    frame = getframe(gcf);
    writeVideo(v,frame);


for tt = 1:25
    for i = 1:numOfP
        uInf{i}.ug = kappa*[-imag(s{i}.x);0*s{i}.x];
    end
    s = squirmer_solver(s, beta ,dt, uInf, qtype,amp,numOfP);
    
    % figure(),
    plot(real(s{1}.x),imag(s{1}.x))
    hold on, plot(real(s{2}.x),imag(s{2}.x))
    plot(real(s{3}.x),imag(s{3}.x))
    hold off
    axis equal
    xlim([-5 20])
    ylim([-8 8])
    pause(0.05)
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);
