function test_random()
%
% function test_random();
%
% Test gkypsdp_solver and compare it to SeDuMi and SDPT3
% 
% Generate random generalized KYP SDP, primal is strictly feasible, dual
% may be unbounded below, so need to run multiple instances to obtain a
% primal / dual feasible test case

clear all;

% generate random generalized KYP SDP
% primal feasible, so maybe unbounded below, need to run many time to get a
% primal / dual feasible SDP, 

prob.p = 2;
prob.L = 2;
prob.n = {15,15};
prob.m = {1,1};
for i=1:prob.L
    prob.nm{i}=prob.n{i}+prob.m{i};
end
Phicase = {1,2}; % 1 = continuous time, 2 = discrete time

L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;

prob.w=randn(p,1);
for i=1:L
    k=0;
    while sum(k)<n{i}
        if Phicase{i}==1
            sys=rss(n{i},m{i},m{i});
        elseif Phicase{i}==2
            sys=drss(n{i},m{i},m{i});
        end
        if m{i}==2
            sys.B(:,2)=zeros(n{i},1);
        end
        %[Abar,Bbar,Cbar,T,k] = ctrbf(sys.A,sys.B,sys.C); 1
        k=rank(ctrb(sys.A,sys.B));
    end
    prob.A{i}=sys.A;
    prob.B{i}=sys.B;
    
    if Phicase{i}==1
        prob.Phi{i}=[0,1;1,0]; 
        if 1==1
            prob.Psi{i}=[0,0;0,0];
        else
            wL = (rand()-0.5)*200;
            wH = (rand()-0.5)*200; %wH = -wL;
            wC = (wL+wH)/2;
            prob.Psi{i} = [-1,j*wC;-j*wC,-wL*wH];
        end
    elseif Phicase{i}==2
        prob.Phi{i}=[1,0;0,-1];
        if 1==1
            prob.Psi{i}=[0,0;0,0];
        else
            wL = (rand()-0.5)*2*pi;
            wH = (rand()-0.5)*2*pi; %wH = -wL;
            wC = (wL+wH)/2;     
            wD = (wH-wL)/2;
            prob.Psi{i} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
        end
    end
    
    % choose x0 and P0 and Q0 random
    x0=randn(p,1);
    P0=randn(n{i},n{i});  P0=P0+P0';
    Q0=randn(n{i},n{i});  [V,D]=eig(Q0*Q0');  Q0=V*diag(abs(randn(n{i},1)))*V';

    % choose M random
    prob.M{i} = zeros(nm{i}^2,p);
    for ii=1:p
       X = randn(nm{i},nm{i});   X = X+X';  
       prob.M{i}(:,ii) = X(:);
    end;

    % subtract random p.d. slack to get N
    S = randn(nm{i},nm{i});   [V,D] = eig(S*S');  S = V*diag(rand(nm{i},1))*V';
    N=-[prob.A{i},prob.B{i};eye(n{i}),zeros(n{i},m{i})]'*(kron(prob.Phi{i},P0)+kron(prob.Psi{i},Q0))...
        *[prob.A{i},prob.B{i};eye(n{i}),zeros(n{i},m{i})]-reshape(prob.M{i}*x0,nm{i},nm{i})-S; 
    %prob.N{i} = 0.5 * (N + N');
    prob.N{i}=eye(size(N))*min(eig(N));
    prob.N{i}(end,end)=-0.1;
end

% use gkypsdp solver
opt.sample = 1;
opt.IPMsolver = 1; 
[sol_1,info_1] = gkypsdp_solver(prob,opt);
[sol_1.Pobj, sol_1.Dobj]
[info_1.iters, info_1.stime/info_1.iters, info_1.ptime]

% solve original sdp using SeDuMi
opt.IPMsolver = 2;
[sol_2,info_2] = gkypsdp_solver(prob, opt);
[sol_2.Pobj, sol_2.Dobj]
[info_2.iters, info_2.stime/info_2.iters, info_2.ptime]

% solve original sdp using SDPT3
opt.IPMsolver = 3;
[sol_3,info_3] = gkypsdp_solver(prob, opt);
[sol_3.Pobj, sol_3.Dobj]
[info_3.iters, info_3.stime/info_3.iters, info_3.ptime]









