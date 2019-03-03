function ex_NORM()
%
% function ex_NORM();
%
% A simple example of computing the H_infinity norm of a transfer function
% The purpose is to illustrate how to enter problem in matlab and call
% gkypsdp_solver.
% Refer to Section 3.2 in the gKYPSDP user guide.

disp('******* A simple example of computing the H_infinity norm ********');

% transfer function: A,b,c,d
A = [-6, -7.5, -12.5; 2, 0, 0; 0, 2, 0];
b = [4;0;0];
c = [-1.5, 0.5, -0.6];
d = 1;

% construct gKYP-SDP
prob.w = [1];
%prob.p = 1;
prob.L = 1;
%prob.n = {3};
%prob.m = {1};
prob.A{1} = A;
prob.B{1} = b;
prob.Phi{1} = [0,1;1,0];
prob.Psi{1} = [0,0;0,0];
prob.N{1} = [c'*c, c'*d; d'*c, d'*d];
Mtmp = [zeros(3,4); zeros(1,3),-1];
prob.M{1}(:,1) = Mtmp(:);

% solve gKYP-SDP
[sol,info] = gkypsdp_solver(prob);

% H infinity norm
gamma = sqrt(sol.x);

% plot transfer function magnitude
figure;
n = 200;
w = logspace(-1,2,n);
Hw = zeros(1,n);
for ii=1:n
    Hw(ii) = c*((j*w(ii)*eye(3)-A)\b)+d;
end
loglog(w,abs(Hw),...
    w,gamma*ones(1,n),'r--');
xlabel('Frequency (radian/second)');
ylabel('Magnitude');
axis([0.1,100,0.35,4]);



