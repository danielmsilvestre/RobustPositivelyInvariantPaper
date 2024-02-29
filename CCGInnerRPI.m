function RPI = CCGInnerRPI(A, U, N)
%CCGINNERRPI Summary of this function goes here
%   Detailed explanation goes here

% We need to get a box set since currently it only supports lp balls
if ~isempty(U.A) || norm(U.G-diag(diag(U.G))) >= 1E-3
    U = boxCCG(U);
end

% value equivalent to infinity
Kinf = 1000;

% Numerically calculate overbound for the remainder of the expression (NEED
% TO BE IMPROVED)
minSVD = 0;
rest = zeros(length(A));
for i = N+1:Kinf
    rest = rest + A^i;
    minSVD = minSVD + min(svd(A^i));
end

n = size(A,1);

%Rest.G = min(svd(rest))*eye(2)*min(diag(U.G)); 
Rest.G = minSVD*eye(2)*min(diag(U.G)); 
Rest.c = zeros(n,1);
Rest.A = zeros(0,n);
Rest.b = zeros(0,1);
Rest.type = 2;
Rest.idx = n;

% Explicitly calculate the RPI set
RPI = U;
for i = 1:N
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,U,zeros(n,1)));
end

RPI = CCGMinkowskiSum(RPI, Rest);

end

