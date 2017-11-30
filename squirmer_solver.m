function [s, mu]= squirmer_solver(s ,beta, dt, uInf,qtype,amp,numOfP)
% solve for rigid body translational velocity and angular velocity
% input 's' is a cell array of structures
% 'uInf' is a cell array of structures with attributes:
%       .ug, background velocity
%       .pg, background pressure
%       'numOfQ' is total number of quadrature points
%       'numOfP' is number of squirmers

%% number of squirmers and quadrature points
% numOfP = numel(s);
numOfQ = 0;
for i = 1:numOfP
    numOfQ = numOfQ + s{i}.n;
end

%% assemble matrix A for numOfP squirmers
As = zeros(2*numOfQ);
incr_i = 0; incr_j = 0;
for i = 1:numOfP
    incr_j = 0;
    for j = 1:numOfP
        if (i==j)
            if qtype=='g'
                Asisj = SLPselfmatrix(s{i});
            else
                [Asisj,~] = stokesselfevalm(s{i}, s{j}, s{j}.n, 's', 'e', 'C');
            end
        else
            Asisj = SLPmatrix2(s{i},s{j});
        end
        As(incr_i+1:incr_i+2*s{i}.n, incr_j+1:incr_j+2*s{j}.n) = Asisj;
        incr_j = incr_j + 2*s{j}.n;
    end
    incr_i = incr_i + 2*s{i}.n;
end

%% assemble matrix B (rigid body velocity)
B = zeros(2*numOfQ,3*numOfP);
incr_i = 0; 
for i = 1:numOfP
    Bi = [-ones(s{i}.n,1),zeros(s{i}.n,1),imag(s{i}.x)-imag(s{i}.xc);...
        zeros(s{i}.n,1),-ones(s{i}.n,1),-(real(s{i}.x)-real(s{i}.xc))];
    B(incr_i+1:incr_i+2*s{i}.n,3*(i-1)+1:3*i) = Bi;
    incr_i = incr_i + 2*s{i}.n;
end

%% assemble matrix C (force and torque free)
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

%% lower right submatrix
Q = zeros(3*numOfP,3*numOfP);
A = [As,B;C,Q];

%% assemble rhs swimming pattern (all same type of squirmer depending on beta)
rhs = zeros(2*numOfQ+3*numOfP,1);
incr_i = 0;
for i = 1:numOfP
    rhs(incr_i+1:incr_i+2*s{i}.n) = amp*[sin(s{i}.t)+beta*sin(2*s{i}.t);...
        sin(s{i}.t)+beta*sin(2*s{i}.t)].*[real(s{i}.tang);imag(s{i}.tang)]...
        -uInf{i}.ug;    % modify right hand side based on background velocity
    incr_i = incr_i + 2*s{i}.n;
end

%% solve for rigid body velocity
soln = (C*(As\B))\(C*(As\rhs(1:end-3*numOfP)));
mu = As\(rhs(1:end-3*numOfP)-B*soln);


%% update position and evolve one step in time
for i = 1:numOfP
    Ux = soln(3*i-2); Uy = soln(3*i-1); w = soln(3*i);
    s{i}.theta = s{i}.theta + w*dt;
    s{i}.xc = s{i}.xc + Ux*dt + 1i*Uy*dt;
    s{i}.Z = @(t) s{i}.xc+cos(s{i}.theta)*s{i}.Zr(t)-sin(s{i}.theta)*s{i}.Zi(t) +...
        1i*( sin(s{i}.theta)*s{i}.Zr(t)+cos(s{i}.theta)*s{i}.Zi(t) );
    [s{i}, N, np] = quadr_pan(s{i}, s{i}.n, 'p', 'C');
    s{i}.n = N;
end
