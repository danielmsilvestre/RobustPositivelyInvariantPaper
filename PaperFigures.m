close all;
clearvars;
clc;

% Check if we need to create a folder for figures
yourFolder = 'PaperFigs';
if ~exist(yourFolder, 'dir')
    mkdir(yourFolder)
end

rng('default'); % comment this to generate other examples


%% Fig.1 - Easy example for a square as initial condition
% value equivalent to infinity
Kinf = 1000;

% horizon H
H = 5;

V = rand(2);
A = V*diag([0.8 0.9])*V^-1;

% Define set for the actuation as a box of size 2
W.G = 2*eye(2);
W.c = zeros(2,1);
W.A = zeros(0,2);
W.b = zeros(0,1);
W.type = inf;
W.idx = 2;

% Explicitly calculate the RPI set
RPI = W;
for i = 1:Kinf
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,W,zeros(2,1)));
end

% RPI set from the iterative method
maxError = 1E-6;
RPIiter = W;
i = 1;
while norm(A^i)*sqrt(2)*max(abs(diag(W.G))) >= maxError
    RPIiter = CCGMinkowskiSum(RPIiter, CCGLinMap(A^i,W,zeros(2,1)));
    i = i + 1;
end

RPIInner = CCGInnerRPI(A,W,H);
RPIOuter = CCGOuterRPI(A,W,H);
RPIOuterApprox = CCGOuterRPIApprox(A,W,H);
RPItight = CCGInnerRPITight(A,W,H);

% Draw All sets
[Frpi,prpi] = compileCCG(RPI);
[FrpiOut,prpiOut] = compileCCG(RPIOuter);
[FrpiOuta,prpiOuta] = compileCCG(RPIOuterApprox);
[FrpiIn,prpiIn] = compileCCG(RPIInner);
[FrpiInt,prpiInt] = compileCCG(RPItight);
h = figure;
plot(FrpiOut,prpiOut,'k');
hold on;
plot(FrpiOuta,prpiOuta,'y');
plot(Frpi,prpi,'r');
plot(FrpiInt,prpiInt,'m');
plot(FrpiIn,prpiIn,'g');

xlabel('$x$-coordinate','Interpreter','latex');
ylabel('$y$-coordinate','Interpreter','latex');
legend('$\Gamma_5^\mathrm{outer}$','$\tilde{F}_\infty$', '$F_\infty$','$\Gamma_5^\mathrm{inner}$','$\Gamma_5^\ell$','Interpreter','latex');

saveas(h, strcat(yourFolder,'/Fig1-rpiH',num2str(H)),'pdf');

fprintf("Results from Fig.1:\n " + ...
    "RPI following iterative method with maxError:%e \n " + ...
    "Sizes:\n\t G [%d x %d], c [%d x 1], A [%d x %d], b [%d x 1].\n",maxError,size(RPIiter.G),size(RPIiter.c,1),size(RPIiter.A),size(RPIiter.b,1));


%% Fig.2 - Check evolution of the Outer RPI for various H

% value equivalent to infinity
Kinf = 1000;

V = rand(2);
A = V*diag([0.8 0.9])*V^-1;

% Define set for the actuation as a box of size 2
W.G = 2*eye(2);
W.c = zeros(2,1);
W.A = zeros(0,2);
W.b = zeros(0,1);
W.type = inf;
W.idx = 2;

% Explicitly calculate the RPI set
RPI = W;
for i = 1:Kinf
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,W,zeros(2,1)));
end


H = [5,10,12,40];
minH = min(H);
maxH = max(H);
nbHorizons = length(H);
RPIOuter = cell(1,nbHorizons);
F = cell(1,nbHorizons);
p = cell(1,nbHorizons);

h = figure;
hold on;
xlabel('$x$-coordinate','Interpreter','latex');
ylabel('$y$-coordinate','Interpreter','latex');

legendArray = cell(1,nbHorizons+1);
colorGradient = 1-logspace(-1,0,nbHorizons);

for i = 1:nbHorizons
    % Calculate the various OuterRPI
    RPIOuter{i} = CCGOuterRPI(A,W,H(i));

    % Draw All sets
    [F{i},p{i}] = compileCCG(RPIOuter{i});
    plot(F{i},p{i},colorGradient(i)*ones(1,3));
    legendArray{i} = strcat('$\Gamma_{',num2str(H(i)),'}^\mathrm{outer}$');
end
legendArray{end} = '$F_\infty$';
[Frpi,prpi] = compileCCG(RPI);
plot(Frpi,prpi,'r');
legend(legendArray{:},'Interpreter','latex');


saveas(h, strcat(yourFolder,'/Fig2-rpiOuterH',num2str(minH),'-',num2str(maxH)),'pdf');


%% Fig.3 - Check evolution of the Inner RPI for various H

% value equivalent to infinity
Kinf = 1000;

V = rand(2);
A = V*diag([0.8 0.9])*V^-1;

% Define set for the actuation as a box of size 2
W.G = 2*eye(2);
W.c = zeros(2,1);
W.A = zeros(0,2);
W.b = zeros(0,1);
W.type = inf;
W.idx = 2;

% Explicitly calculate the RPI set
RPI = W;
for i = 1:Kinf
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,W,zeros(2,1)));
end


H = [12,6,4,1];
minH = min(H);
maxH = max(H);
nbHorizons = length(H);
RPIInner = cell(1,nbHorizons);
F = cell(1,nbHorizons);
p = cell(1,nbHorizons);

h = figure;
hold on;
xlabel('$x$-coordinate','Interpreter','latex');
ylabel('$y$-coordinate','Interpreter','latex');

legendArray = cell(1,nbHorizons+1);
colorGradient = logspace(-1,0,nbHorizons);

legendArray{1} = '$F_\infty$';
[Frpi,prpi] = compileCCG(RPI);
plot(Frpi,prpi,'r');


for i = 1:nbHorizons
    % Calculate the various InnerRPI
    RPIInner{i} = CCGInnerRPITight(A,W,H(i));

    % Draw All sets
    [F{i},p{i}] = compileCCG(RPIInner{i});
    plot(F{i},p{i},colorGradient(i)*ones(1,3));
    legendArray{i+1} = strcat('$\Gamma_{',num2str(H(i)),'}^\mathrm{inner}$');
end
legend(legendArray{:},'Interpreter','latex');


saveas(h, strcat(yourFolder,'/Fig3-rpiInnerH',num2str(minH),'-',num2str(maxH)),'pdf');

%% Fig.4 - Check evolution of the Inner RPI for various H and random W

% value equivalent to infinity
Kinf = 500;

V = rand(2);
A = V*diag([0.8 0.9])*V^-1;
B = eye(2);

% Define the random set for the actuation
W.G = randn(2,20);
W.c = randn(2,1);
W.A = randn(10,20);
W.b = randn(10,1);
W.type = [inf,2];
W.idx = 10*ones(1,2);

% Explicitly calculate the RPI set
fprintf("Calculating a very high-precision RPI for drawing purposes (really slow)...\n");
RPI = W;
for i = 1:Kinf
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,W,zeros(2,1)));
end

% RPI set from the iterative method
maxError = 1E-6;
fprintf("Calculating RPI through iterative method (may take a while)...\n");
RPIiter = W;
i = 1;
tic;
while norm(A^i)*sqrt(2)*max(abs(sum(W.G,2))) >= maxError
    RPIiter = CCGMinkowskiSum(RPIiter, CCGLinMap(A^i,W,zeros(2,1)));
    i = i + 1;
end
computingTimeIterativeMethod = toc;


H = [12,6,4,1];
minH = min(H);
maxH = max(H);
nbHorizons = length(H);
RPIInner = cell(1,nbHorizons);
F = cell(1,nbHorizons);
p = cell(1,nbHorizons);

h = figure;
hold on;
xlabel('$x$-coordinate','Interpreter','latex');
ylabel('$y$-coordinate','Interpreter','latex');

legendArray = cell(1,nbHorizons+1);
colorGradient = logspace(-1,0,nbHorizons);

legendArray{1} = '$F_\infty$';
[Frpi,prpi] = compileCCG(RPI);
plot(Frpi,prpi,'r');

computingTimeH = zeros(1,nbHorizons);

for i = 1:nbHorizons
    % Calculate the various InnerRPI
    tic;
    RPIInner{i} = CCGInnerRPITight(A,W,H(i));
    computingTimeH(i) = toc;

    % Draw All sets
    [F{i},p{i}] = compileCCG(RPIInner{i});
    plot(F{i},p{i},colorGradient(i)*ones(1,3));
    legendArray{i+1} = strcat('$\Gamma_{',num2str(H(i)),'}^\mathrm{inner}$');
end
legend(legendArray{:},'Interpreter','latex');


saveas(h, strcat(yourFolder,'/Fig4-rpiInnerH',num2str(minH),'-',num2str(maxH)),'pdf');
fprintf("Results from Fig.4:\n RPI following iterative method with maxError:%e \n " + ...
    "Sizes:\n\t G [%d x %d], c [%d x 1], A [%d x %d], b [%d x 1].\n" + ...
    "Computing Time for Iterative Method: %e seconds \n" + ...
    "Computing time for each horizon: [ %g ]",maxError,size(RPIiter.G),size(RPIiter.c,1),size(RPIiter.A),size(RPIiter.b,1), computingTimeIterativeMethod, computingTimeH);

fprintf("Sizes of the Inner RPIs for each horizon H:\n")
for i = 1:nbHorizons
    fprintf("Horizon: %d \t G [%d x %d], c [%d x 1], A [%d x %d], b [%d x 1].\n", H(i) ,size(RPIInner{i}.G),size(RPIInner{i}.c,1),size(RPIInner{i}.A),size(RPIInner{i}.b,1));
end


%% Fig.5 and Fig.6 - Computing time and Volume ratio between RPI and Inner Approximation RPI for 100 systems

nbSystems = 100;
maxH = 20;

volRPI = -ones(1,nbSystems);
volInnerRPI = -ones(nbSystems,maxH);
volComputingTimeRPI = -ones(1,nbSystems);
volComputingTimeInnerRPI = -ones(nbSystems,maxH);

for j = 1:nbSystems
    % Generate random system in ss format
    A = eye(2);
    while norm(A) >= 0.95
        sys = drss(2);
        [A,~,~,~] = ssdata(sys);
    end

    % Define the random set for the actuation
    W = [];
    while isEmptyW(W)
        W.G = randn(2,20);
        W.c = randn(2,1);
        W.A = randn(10,20);
        W.b = randn(10,1);
        W.type = [inf,2];
        W.idx = 10*ones(1,2);
    end
    if mod(j,10) == 0 || j == 1
        fprintf("Calculating RPI through iterative method for %d case (may take a while)...\n",j);
    end
    RPIiter = W;
    i = 1;
    while norm(A^i)*sqrt(2)*max(abs(sum(W.G,2))) >= maxError
        RPIiter = CCGMinkowskiSum(RPIiter, CCGLinMap(A^i,W,zeros(2,1)));
        i = i + 1;
    end
    % Compute volume of the RPI using iterative method
    tic;
    volRPI(j) = CCGVolume(RPIiter);
    volComputingTimeRPI(j) = toc;

    % Check the volume for the various inner sets depending on H
    for H = 1:maxH

        % Calculate the inner RPI set
        innerRPI = CCGInnerRPITight(A,W,H);

        % Compute the volume of the inner RPI
        tic;
        volInnerRPI(j,H) = CCGVolume(innerRPI);
        volComputingTimeInnerRPI(j,H) = toc;

    end
end

h = figure;
hold on;
xlabel('$H$','Interpreter','latex');
ylabel('$\frac{\mathrm{Vol}(\Gamma_H^\mathrm{inner})}{\mathrm{Vol}(F_\infty)}$','Interpreter','latex');

volRatios = volInnerRPI./volRPI';

plot(1:maxH,volRatios,'Color',0.7*ones(1,3));
plot(1:maxH,min(volRatios),'k','LineWidth',2);
plot(1:maxH,max(volRatios),'k','LineWidth',2);
saveas(h, strcat(yourFolder,'/Fig5-volH',num2str(minH),'-',num2str(maxH)),'pdf');

h = figure;
hold on;
xlabel('$H$','Interpreter','latex');
ylabel('elapsed time [s]','Interpreter','latex');
ylim([0,10])

plot(1:maxH,mean(volComputingTimeInnerRPI),1:maxH,mean(volComputingTimeRPI)*ones(1,maxH),'LineWidth',2);
legend('$\Gamma_H^\mathrm{inner}$','$F_\infty$','Interpreter','latex');
saveas(h, strcat(yourFolder,'/Fig6-ctimeH',num2str(minH),'-',num2str(maxH)),'pdf');


%%% Auxiliary Functions %%%

% Function used to test if the random W is empty - it is not meant to be a
% generic test since if we suspect a set to be empty we can simply generate
% a new random set. Such a solution is more efficient for the purpose of
% this script.

function bool = isEmptyW(W)
    if isempty(W)
        bool = true;
        return;
    end
    
    % compute the minimum norm solution of A*xi = b
    xi_sol = (W.A'*((W.A*W.A')\W.b));

    % check if it satisfies the norm bounds for xi
    notEmpty = true;
    i = 1;
    indices = [1, cumsum(W.idx)];
    while notEmpty && i <= length(W.type)
        notEmpty = notEmpty & (norm(xi_sol(indices(i):indices(i+1),1),W.type(i)) <= 1);
        i = i + 1;
    end
    bool = ~notEmpty;
end


