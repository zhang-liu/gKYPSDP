function [sol,info] = origkypsdp(prob,opt)
% ORIGKYPSDP Solve the generalized KYP SDP of the original inequality form
% without problem reformulation
%    minimize     w'*x
%    subject to   [A{i},B{i};I,0]'*(kron(Phi{i},P(i})+kron(Psi{i},Q{i}))*
%                 [A{i},B{i};I,0] + M(x) + N <= 0, i=1..L.
% 
% ORIGKYPSDP(PROB), for problem data PROB, is the optimal solution.
% ORIGKYPSDP(PROB,OPT), for problem data PROB and options OPT, is the
%   optimal solution.

% Problem data:
% - prob.w           px1 vector
% - prob.A{i}        n{i} x n{i} matrix
% - prob.B{i}        n{i} x m{i} matrix
% - prob.Phi{i}      2x2 hermitian matrix
% - prob.Psi{i}      2x2 hermitian matrix
% - prob.M{i,j}      (n+m)x(n+m) hermitian matrix, j=1..p
% - prob.N{i}        (n+m)x(n+m) hermitian matrix
% - prob.p           number of elements in x
% - prob.L           number of constraints
% - prob.n{i}        number of state in i
% - prob.m{i}        number of input in i
% - prob.nm{i}       n{i}+m{i}
%
% Options:
% - opt.IPMsolver   use SeDuMi if 2, use SDPT3 if 3
%
% Solutions: 
% - sol.x          px1 primal vector
% - sol.Pobj       primal objective, w'*x
% - sol.Dobj       dual objective
%
% Solver information:
% - info.ptime      preprocessing CPU time in seconds
% - info.stime      IPM solving CPU time in seconds
% - info.iters      number of iterations

prob.w = real(prob.w);
for i=1:prob.L
    prob.A{i} = real(prob.A{i});
    prob.B{i} = real(prob.B{i});
    prob.N{i} = 0.5*(prob.N{i}+prob.N{i}');
    prob.nm{i} = prob.n{i}+prob.m{i};
    for jj = 1:prob.p
        Mt = 0.5*(reshape(prob.M{i}(:,jj),prob.nm{i},prob.nm{i})+reshape(prob.M{i}(:,jj),prob.nm{i},prob.nm{i})');
        prob.M{i}(:,jj)=Mt(:);
    end
end

L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;
Phi=prob.Phi;

x = sdpvar(p,1);
F = [];
for i=1:L;
P{i} = sdpvar(n{i},n{i},'hermitian','complex');
Q{i} = sdpvar(n{i},n{i},'hermitian','complex');

F = F+set([prob.A{i},prob.B{i};eye(n{i}),zeros(n{i},m{i})]'*(kron(prob.Phi{i},P{i})...
    +kron(prob.Psi{i},Q{i}))*[prob.A{i},prob.B{i};eye(n{i}),zeros(n{i},m{i})]...
    +reshape(prob.M{i}*x,nm{i},nm{i}) + prob.N{i} <= 0) + set(Q{i} >= 0);
end

obj = prob.w'*x;
if opt.IPMsolver == 2
    status = solvesdp(F,obj,sdpsettings('solver','sedumi','savesolveroutput',1));

    sol.x = double(x);
    sol.Dobj=0;
    for i=1:L
        sol.Dobj=sol.Dobj+trace(prob.N{i}*double(dual(F(2*i-1))));
    end
    sol.Pobj = prob.w'*sol.x; 
    info.iters = status.solveroutput.info.iter;
    info.ptime = status.solveroutput.info.timing(1);
    info.stime = status.solveroutput.info.timing(2);
elseif opt.IPMsolver == 3
    status = solvesdp(F,obj,sdpsettings('solver','sdpt3','savesolveroutput',1));
   
    sol.x = double(x);
    sol.Dobj = 0;
    for i=1:L
        sol.Dobj=sol.Dobj+trace(prob.N{i}*double(dual(F(2*i-1))));
    end
    sol.Pobj = prob.w'*sol.x;
    info.iters = status.solveroutput.info.iter;
    info.ptime = status.solveroutput.info.ttime.preproc;
    info.stime = status.solveroutput.info.ttime.pred;
end
