function RPI = CCGInnerRPITight(A, U, N)
%CCGINNERRPI Summary of this function goes here
%   Detailed explanation goes here

% We need to get a box set since currently it only supports lp balls
% if ~isempty(U.A) || norm(U.G-diag(diag(U.G))) >= 1E-3
%     U = boxCCG(U);
% end

rest = inv(eye(size(A))-A);
for i = 0:N
    rest = rest-A^i;
end

n = size(A,1);

Rest = CCGLinMap(rest,U,zeros(n,1));
% Rest.G = rest*U.G;
% Rest.c = zeros(n,1);
% Rest.A = zeros(0,n);
% Rest.b = zeros(0,1);
% Rest.type = U.type;
% Rest.idx = n;

% Explicitly calculate the RPI set
RPI = U;
for i = 1:N
    RPI = CCGMinkowskiSum(RPI, CCGLinMap(A^i,U,zeros(n,1)));
end

RPI = CCGMinkowskiSum(RPI, Rest);

end

