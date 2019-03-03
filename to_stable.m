function [prob] = to_stable(prob,rho)
%
% function [prob] = to_stable(prob,rho);
%
% Perform a LQR state feedback to transform KYP LMI so that A{i} is stable

L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;
Phi=prob.Phi;

for i=1:L
    if n{i}~=0
        eigA=eig(prob.A{i});
        if max(real(Phi{i}(1,1)*eigA.*conj(eigA)+Phi{i}(1,2)*conj(eigA)+Phi{i}(2,1)*eigA+Phi{i}(2,2)))>= - 1E-1
            Q=eye(n{i});
            R=rho*eye(m{i});
            if sum(sum(abs(Phi{i}-[0,1;1,0])))<1E-6
                [F,S,e]=lqr(prob.A{i},prob.B{i},Q,R);
            elseif sum(sum(abs(Phi{i}-[1,0;0,-1])))<1E-6
                [F,S,e]=dlqr(prob.A{i},prob.B{i},Q,R);
            end
            prob.A{i}=prob.A{i}-prob.B{i}*F;
            prob.Kt{i}=prob.Kt{i}*[eye(n{i}),zeros(n{i},m{i});-F,eye(m{i})];
        end
    end
end
