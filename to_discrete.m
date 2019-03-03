function [prob] = to_discrete(prob)
%
% function [prob] = to_discrete(prob)
%
% Transform generalized KYP LMI into the discrete-time KYP LMI

L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;

for i=1:L
    if n{i}~=0
        [U,T]=schur(prob.Phi{i});
        if T(1,1) < 0
            U=U*[0,1;1,0];
            T=[0,1;1,0]*T*[0,1;1,0];
        end
        U=[sqrt(T(1,1)),0;0,sqrt(-T(2,2))]*U';
        prob.Phi{i}=[1,0;0,-1];  % Phi = U'*T*U;
        prob.Psi{i}=inv(U')*prob.Psi{i}*inv(U);
        Tt=inv(U(2,1)*prob.A{i}+U(2,2)*eye(n{i}));
        prob.Kt{i}=prob.Kt{i}*[Tt,-Tt*U(2,1)*prob.B{i};zeros(m{i},n{i}),eye(m{i})];
        prob.A{i}=(U(1,1)*prob.A{i}+U(1,2)*eye(n{i}))*Tt;
        prob.B{i}=U(1,1)*prob.B{i}-prob.A{i}*(U(2,1)*prob.B{i});
    end
end
