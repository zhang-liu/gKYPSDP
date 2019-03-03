function [sol,info] = gkypsdp_solver(prob,opt)
% GKYPSDP_SOLVER Solve the generalized KYP SDP of the form
%    minimize     w'*x
%    subject to   [A{i},B{i};I,0]'*(kron(Phi{i},P(i})+kron(Psi{i},Q{i}))*
%                 [A{i},B{i};I,0] + M(x) + N <= 0, i=1..L.
%
% GKYPSDP_SOLVER(PROB), for problem data PROB, is the optimal solution.
% GKYPSDP_SOLVER(PROB,OPT), for problem data PROB and options OPT, is the
%   optimal solution.

% Problem data:
% - prob.w              px1 vector
% - prob.A{i}           n{i} x n{i} matrix
% - prob.B{i}           n{i} x m{i} matrix
% - prob.Phi{i}         2x2 hermitian matrix
% - prob.Psi{i}         2x2 hermitian matrix
% - prob.M{i,j}         (n+m)x(n+m) hermitian matrix (see user guide)
% - prob.N{i}           (n+m)x(n+m) hermitian matrix
% - prob.p              number of elements in x
% - prob.L              number of constraints
% - prob.n{i}           number of state in i
% - prob.m{i}           number of input in i
%
% Options:
% - opt.sample      sampling method, default = 1
% - opt.IPMsolver   interior-point method solver, default = 1
% - opt.feasx       initial value of x, default = 0
% - opt.roh         roh for LQR feedback, default = 0.25
% - opt.maxiters    number of iterations, default = 50
% - opt.abstol      absolution tolerance of duality gap, default = 1E-6
% - opt.reltol      relative tolerance of duality gap, default = 1E-6
% - opt.pfeastol    residual tolerance in the primal feasibility equations, default = 1E-6
% - opt.dfeastol    residual tolerance in the dual feasibility equations, default = 1E-6
%
% Solutions:
% - sol.x               px1 primal variable at last iteration
% - sol.z               reformulated dual variable at last iteration
% - sol.Pobj            primal objective, w'*x
% - sol.Dobj            dual objective, 
%
% Solver information:
% - info.ptime          preprocessing CPU time in seconds
% - info.stime          IPM solving CPU time in seconds
% - info.iters          number of iterations
% - info.pfeas          residual in the primal feasibility equations
% - info.dfeas          residual in the dual feasibility equations
% - info.absgap         absolution gap
% - info.relgap         relative gap 
%
% Please send any feedback to zhang@ee.ucla.edu | vandenbe@ee.ucla.edu
 
ptimestart = cputime;

prob.p = size(prob.w,1);
for i=1:prob.L
    prob.n{i} = size(prob.B{i},1);
    prob.m{i} = size(prob.B{i},2);
end

% solve the original SDP using SeDuMi or SDPT3
if nargin == 2 & isfield(opt, 'IPMsolver')==1 & opt.IPMsolver > 1
    [sol,info] = origkypsdp(prob, opt);
    return;
end

% process the problem data
prob.w = real(prob.w);
for i=1:prob.L
    prob.nm{i} = prob.n{i}+prob.m{i};
    if size(prob.A{i},1)~=prob.n{i} | size(prob.A{i},2)~=prob.n{i}
        error('size of prob.A{i} does not match prob.n{i}');
    end
    prob.A{i} = real(prob.A{i});
    if size(prob.B{i},1)~=prob.n{i} | size(prob.B{i},2)~=prob.m{i}
        error('size of prob.B{i} does not match prob.n{i} or prob.m{i}');
    end
    prob.B{i} = real(prob.B{i});
    if size(prob.N{i},1)~=prob.nm{i} | size(prob.N{i},2)~=prob.nm{i}
        error('size of prob.N{i} does not match prob.n{i}, prob.m{i}');
    end
    prob.N{i} = sparse(0.5*(prob.N{i}+prob.N{i}'));
    if size(prob.M{i},1)~=(prob.nm{i})^2 | size(prob.M{i},2)~=prob.p
        error('size of prob.M{i} does not match prob.mn{i}, prob.p{i}');
    end
    for jj = 1:prob.p
        Mt = 0.5*(reshape(prob.M{i}(:,jj),prob.nm{i},prob.nm{i})+reshape(prob.M{i}(:,jj),prob.nm{i},prob.nm{i})');
        prob.M{i}(:,jj)=sparse(Mt(:));
    end
    prob.Kt{i} = eye(prob.nm{i},prob.nm{i});
end

% this is for convenience
L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;
Phi=prob.Phi;

% set default value for options if it is not specified
if nargin < 2 | isfield(opt, 'sample')==0
    opt.sample = 1;
end
if nargin < 2 | isfield(opt, 'IPMsolver')==0
    opt.IPMsolver = 1;
end
if nargin < 2 | isfield(opt, 'roh')==0
    opt.roh = 0.25;
end
if nargin < 2 | isfield(opt, 'maxiters')==0
    opt.maxiters = 50;
end
if nargin < 2 | isfield(opt, 'abstol')==0
    opt.abstol = 1E-6;
end
if nargin < 2 | isfield(opt, 'reltol')==0
    opt.reltol = 1E-6;
end
if nargin < 2 | isfield(opt, 'pfeastol')==0
    opt.pfeastol = 1E-6;
end
if nargin < 2 | isfield(opt, 'dfeastol')==0
    opt.dfeastol = 1E-6;
end
if nargin < 2 | isfield(opt, 'feasx')==0
    % heuristic for determining initial x, not necessary feasible
    ii=1;
    MM = [];
    for i=1:L
        if m{i}==2
            MM(ii,:) = prob.M{i}(end,:);
            NN(ii,:) = prob.N{i}(end,end);
            if NN(ii,:) < 0
                none(ii,:) = NN(ii,:);
            else
                none(ii,:) = -0.5; 
            end
            ii=ii+1;
        end
    end
    opt.feasx = zeros(p,1);
    if isempty(MM)~=1 & rank(MM) > 0
        opt.feasx = MM\(none-NN);
    end
end
   
% check if the problem is within the solver limitation
for i=1:L
    if n{i}==0
        error('n{i}==0, currently only support n{i}>0'); 
    elseif m{i}>2 | (m{i}==2 & sum(abs(prob.B{i}(:,2)))>0)
        error('Currently only support single input B=b or B=[b,0]'); 
    elseif sum(sum(abs(Phi{i}-[1,0;0,-1])))>1E-6 & sum(sum(abs(Phi{i}-[0,1;1,0])))>1E-6
        error('Currently only support continuous or discrete Phi'); 
    elseif rank(ctrb(prob.A{i},prob.B{i}))<n{i}
        error('(A{i},B{i}) is not controllable'); 
    end
end
      
% solve the gKYPSDP
[prob] = to_stable(prob,opt.roh);    
[prob] = to_discrete(prob);
[sos_data,err] = to_sos(prob, opt.sample,opt.feasx);
%[prob.F2] = to_stable(prob,0.25);
if err==1
    sol = 0; info.err=1; return; 
end
info.ptime = cputime - ptimestart;
[y,X,z,ptime,stime,iters,p_feas,d_feas,abs_gap,rel_gap]=solSOSSDP(...
    sos_data.As,sos_data.Bs,sos_data.bs,sos_data.Cs,sos_data.qs,sos_data.Qs,opt,sos_data.varargin);
 
% process output
sol.x = y;
sol.z = z;
sol.Pobj = prob.w'*sol.x;
sol.Dobj = 0;
for i=1:prob.L
    sol.Dobj = sol.Dobj + sos_data.bs{i}'*z{i};
end
info.ptime = info.ptime + ptime;
info.stime = stime;
info.iters = iters;
info.pfeas = p_feas;
info.dfeas = d_feas;
info.absgap = abs_gap;
info.relgap = rel_gap;

