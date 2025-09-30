% Example from Section 6.2 of Vinod, Weiss, Di Cairano, "Projection-free 
% computation of robust controllable sets with constrained zonotopes", 
% Automatica, 2025.

% Computing the T-step Robust Controllable (RC) set
T = 20; % Number of time steps

A = [0.99 0.02; -0.15 0.99];
B = [-0.01; 0.08];
F = eye(2);

U = zono(1.5,0);
W = zono(0.01*eye(2),zeros(2,1));
Goal = zono(0.5*eye(2),[1.5; 0]);
X = zono(10*eye(2),zeros(2,1));
X = halfspaceIntersection(X,[-1 0;2 1],[2;5]);

figure; hold on
plot(X,'y',1)
xlim([-3 5])
ylim([-10 6])
xlabel('x_1')
ylabel('x_2')

ratiog1=[];
ratioc1=[];
cum_time=[];

InnerApproxFlag = 1;
N = 100;
ReduceFlagEveryNSteps = 1;
ReduceFlagAtEnd = 1;
tStart = tic;
Ainv = inv(A);
K = conZono(Goal);
tcum=tic;
for i = 1:T
    if InnerApproxFlag
        % Approximate from paper above
        W0 = W + (-W.c);
        Gamma = ([K.G;K.A])\[eye(K.n);zeros(K.nC,K.n)];
        D = zeros(K.nG);
        for j = 1:K.nG
            ej = zeros(K.nG,1);
            ej(j) = 1;
            dir = ej'*Gamma;
            D(j,j) = 1 - sum(abs(dir*W0.G));
        end
        KminusWInner = conZono(K.G*D,K.c-W.c,K.A*D,K.b);
        % figure;
        % plot(KminusW,'r',0.1)
        % plot(KminusWInner,'b',0.1)
        % drawnow
        K = Ainv*(KminusWInner + (-B*U));
    else
        KminusW = pontryDiff(K,W);
        K = Ainv*(KminusW + (-B*U));
    end
    
    K = halfspaceIntersection(K,[-1 0;2 1],[2;5]);
    nGpre = K.nG;
    nCpre = K.nC;

    ratiog1(end+1,:)=K.nG;
    ratioc1(end+1,:)=K.nC;

    if ReduceFlagEveryNSteps
        if mod(i,N) == 0
            K = reduceConZonoGurobi(K);          
        end
    end
    cum_time(end+1,:)=toc(tcum);
end
t_Calc = toc(tStart)
tStart = tic;
if ReduceFlagAtEnd
    K = reduceConZonoGurobi(K);
    ratiog1(end+1,:)=K.nG;
    ratioc1(end+1,:)=K.nC;
    cum_time(end+1,:)=toc(tcum);
end
t_Reduce = toc(tStart)
[i nGpre K.nG nCpre K.nC]

tStart = tic;
plot(K,'m',1)
t_Plot = toc(tStart)
plot(Goal,'r',1)
drawnow

%% Initial condition containment check
% X0 = zono(eye(2),[1;-5]);
% [verts,~] = plot(X0,'g',1);

% tStart = tic;
% nVerts = size(verts,2);
% pointContained = zeros(nVerts,1);
% for k = 1:nVerts
%     pointContained(k) = checkPointContain(K,verts(k,:)');
% end
% setContained = min(pointContained);
% t_setContained = toc(tStart)

%% Evaluate multiple support functions
tStart = tic;
nDirs = 500;
infeasCount = 0;
for k = 1:nDirs
    rng(k)
    [s,x] = supportFunc(K,rand(2,1));
    if size(s,2) ~= 1
        infeasCount = infeasCount + 1;
    end
end
t_support = toc(tStart)
infeasCount


%% Could try to see if set has converged
