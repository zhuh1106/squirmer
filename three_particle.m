function thress_particle

clear all
close all

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

N1 = 400;
N2 = 400;
N3 = 300;
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

% rhs condition
numOfP = 3; %num of squirmer
numOfQ = N1+N2+N3; %num of quadrature
rhs = zeros(2*numOfQ+3*numOfP,1);
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
    As = zeros(2*numOfQ);
    incr_i = 0; incr_j = 0;
    for i = 1:numOfP
        incr_j = 0;
        for j = 1:numOfP
            if (i==j)
                [Asisj,~] = stokesselfevalm(s{i}, s{j}, N1, 's', 'e', 'C');
            else
                Asisj = SLPmatrix2(s{i},s{j});
            end
            As(incr_i+1:incr_i+2*s{i}.n, incr_j+1:incr_j+2*s{j}.n) = ...
                    Asisj;
            incr_j = incr_j + 2*s{j}.n;
        end
        incr_i = incr_i + 2*s{i}.n;
    end
    
    B = zeros(2*numOfQ,3*numOfP);
    incr_i = 0; 
    for i = 1:numOfP
        Bi = [-ones(s{i}.n,1),zeros(s{i}.n,1),imag(s{i}.x)-imag(s{i}.xc);...
            zeros(s{i}.n,1),-ones(s{i}.n,1),-(real(s{i}.x)-real(s{i}.xc))];
        B(incr_i+1:incr_i+2*s{i}.n,3*(i-1)+1:3*i) = Bi;
        incr_i = incr_i + 2*s{i}.n;
    end
    
    C = zeros(3*numOfP,2*numOfQ);
    incr_i = 0; 
    for i = 1:numOfP
        Api = SLPmatrixp(s{i},s{i});
        Cfi = [s{i}.ws',zeros(1,s{i}.n); zeros(1,s{i}.n),s{i}.ws']*Api;
        Cti = [-(imag(s{i}.x)-imag(s{i}.xc))',(real(s{i}.x)-real(s{i}.xc))'];
        Ci = [Cfi;Cti];
        C(3*(i-1)+1:3*i,incr_i+1:incr_i+2*s{i}.n) = Ci;
        incr_i = incr_i + 2*s{i}.n;
    end
    
    Q = zeros(3*numOfP,3*numOfP);
    A = [As,B;C,Q];
    incr_i = 0;
    for i = 1:numOfP
        rhs(incr_i+1:incr_i+2*s{i}.n) = 20*[sin(s{i}.t)+beta*sin(2*s{i}.t);...
            sin(s{i}.t)+beta*sin(2*s{i}.t)].*[real(s{i}.tang);imag(s{i}.tang)];
        incr_i = incr_i + 2*s{i}.n;
    end
    
    soln = (C*(As\B))\(C*(As\rhs(1:end-3*numOfP)))
    for i = 1:numOfP
        Ux = soln(3*i-2); Uy = soln(3*i-1); w = soln(3*i);
        s{i}.theta = s{i}.theta + w*dt;
        s{i}.xc = s{i}.xc + Ux*dt + 1i*Uy*dt;
        s{i}.Z = @(t) s{i}.xc+cos(s{i}.theta)*s{i}.Zr(t)-sin(s{i}.theta)*s{i}.Zi(t) +...
                1i*( sin(s{i}.theta)*s{i}.Zr(t)+cos(s{i}.theta)*s{i}.Zi(t) );
        [s{i}, N, np] = quadr_pan(s{i}, s{i}.n, 'p', 'C');
        s{i}.n = N;
    end
    
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
