function squirmer_main
clear all
close all

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

pos = zeros(3,10);

for i = 1:10

N1 = i*100;
N2 = i*100;
N3 = i*100;
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
[s{1}, N1, np] = quadr_pan(s{1}, N1, 'p', 'C');
s{1}.n = N1;
s{1}.theta = 0;

s{2}.Z  = @(t) x2c + 2*cos(t) + 1i*sin(t);
s{2}.Zr = @(t) real(s{2}.Z(t) - x2c);
s{2}.Zi = @(t) imag(s{2}.Z(t) - x2c);
s{2}.xc = x2c;
[s{2}, N2, np] = quadr_pan(s{2}, N2, 'p', 'C');
s{2}.n = N2;
s{2}.theta = 0;

s{3}.Z  = @(t) x3c + 2*cos(t) + 1i*sin(t);
s{3}.Zr = @(t) real(s{3}.Z(t) - x3c);
s{3}.Zi = @(t) imag(s{3}.Z(t) - x3c);
s{3}.xc = x3c;
[s{3}, N3, np] = quadr_pan(s{3}, N3, 'p', 'C');
s{3}.n = N3;
s{3}.theta = 0;

numOfP = 3; %num of squirmer
numOfQ = N1+N2+N3; %num of quadrature
rhs = zeros(2*numOfQ+3*numOfP,1);
beta = -0.4;

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
    s = squirmer_solver(s, beta ,dt);
    % figure()
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
pos(:,i) = [s{1}.xc,s{2}.xc,s{3}.xc];
end
err = pos - pos(:,end)*ones(1,10);
figure(), semilogy(100*(1:9),err(1,1:end-1),'o')
figure(), semilogy(100*(1:9),err(2,1:end-1),'o')
figure(), semilogy(100*(1:9),err(3,1:end-1),'o')
