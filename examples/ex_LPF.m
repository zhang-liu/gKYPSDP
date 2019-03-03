function ex_LPF()
%
% function ex_LPF();
%
% Low pass filter (LPF) design
% User can specify different filter design parameters in ex_LPF.m
% Refer to Section 4.4 in the gKYPSDP user guide.

disp('*************** Low pass filter (LPF) design ********************');

% filter design paramters
n = 50;
wp = 1.0;
ws = 1.3;
ts = 0.01;
d = 10;

% construct gKYP-SDP
prob.w = [zeros(n+1,1);1];
%prob.p = n+2;
prob.L = 2;
%prob.n = {n,n};
%prob.m = {2,2};

II = eye(n);
for i = 1:n
    Mt = sparse([zeros(n+1,n+1),[II(:,i);0]; II(i,:),0,0]);
    Mn(:,i) = Mt(:);
end
Mt = sparse([zeros(n+1,n+1),[zeros(n,1);1];zeros(1,n),1,0]); 
Mn(:,n+1) = Mt(:);

% H(w)^H H(w) <= ts^2,      ws <= w <= pi
wL=ws;   wH=pi;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{1} = [1,0;0,-1];
prob.Psi{1} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{1} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{1} = [[1;zeros(n-1,1)],zeros(n,1)];
Mt = sparse(zeros(n+2,n+2));
prob.M{1} = [Mn Mt(:)];
prob.N{1} = sparse(diag([zeros(n,1);-ts^2;-1]));

% |H(w)-e^{-jdw}|^2 <= tp,  0 <= w <= wp
wL=0;   wH=wp;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{2} = [1,0;0,-1];
prob.Psi{2} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{2} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{2} = [[1;zeros(n-1,1)],zeros(n,1)];
Mt = sparse(diag([zeros(n,1);-1;0]));
prob.M{2} = [Mn Mt(:)];
prob.N{2} = sparse([zeros(n+1,n+1),[-II(:,d);0];-II(d,:),0,-1]);

if 1==1
    % special-purpose solver, gkypsdp
    [sol,info] = gkypsdp_solver(prob);
else
    % general-purpose solver, sedumi
    opt.IPMsolver = 2;
    [sol,info] = gkypsdp_solver(prob,opt);
end

x0 = sol.x(n+1);
x1n = sol.x(1:n);
tp = (sol.x(n+2))^0.5;

wN = 1000;
wr = linspace(0,pi,wN).';
Hwr = repmat(x0,wN,1) + (repmat(exp(-j*wr),1,n).^repmat([1:n],wN,1))*x1n;

% plot filter gain
figure;
semilogy([0,0,wp,wp,0],[1-tp,1+tp,1+tp,1-tp,1-tp],'r--',...
    [ws,ws,pi],[ts*0.1,ts,ts],'r--');
hold on;
semilogy(wr,abs(Hwr));
hold off;
xlabel('w (radian)');
ylabel('Magnitude of H(w)');
axis([0,pi,ts*0.1,2]);

% plot filter phase
figure;
plot([0,2*wp],[0,-d*2*wp],'r--');
hold on;
plot(wr,unwrap(angle(Hwr)));
hold off;
xlabel('w (radian)');
ylabel('Phase of H(w) (radian)');
axis([0,pi,-16,0]);

