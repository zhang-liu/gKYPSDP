function ex_BPF()
%
% function ex_BPF();
%
% Linear-phase band pass filter (BPF) design
% User can specify different filter design parameters in ex_BPF.m 
% Refer to Section 4.3 in the gKYPSDP user guide.

disp('********** Linear-phase band pass filter (BPF) design ************');

% filter design parameters
n=50;
ws1 = 1.6;
ws2 = 2.1;
wp1 = 1.7;
wp2 = 2.0;
tp = 0.05;

% construct gKYP-SDP
prob.w = [zeros(n+1,1);1];
%prob.p = n+2;
prob.L = 5;
%prob.n = {n,n,n,n,n};
%prob.m = {1,1,1,1,1};

II = eye(n);
for i = 1:n
    Mt = sparse([zeros(n,n) II(:,i); II(i,:) 0]);
    Mn(:,i) = Mt(:);
end
Mt = sparse([zeros(n,n+1); zeros(1,n) 1]);
Mn(:,n+1) = Mt(:);

% H(w) >= -ts,      0 <= w <= pi
wL=0;   wH=pi;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{1} = [1,0;0,-1];
prob.Psi{1} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{1} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{1} = [1;zeros(n-1,1)];
Mt = sparse([zeros(n,n+1); zeros(1,n) -1]);
prob.M{1} = [-Mn Mt(:)];
prob.N{1} = sparse(zeros(n+1,n+1));

% H(w) <= ts,       0 <= w <= ws1
wL=0;   wH=ws1;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{2} = [1,0;0,-1];
prob.Psi{2} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{2} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{2} = [1;zeros(n-1,1)];
Mt = sparse([zeros(n,n+1); zeros(1,n) -1]);
prob.M{2} = [Mn Mt(:)];
prob.N{2} = sparse(zeros(n+1,n+1));

% H(w) <= ts,       ws2 <= w <= pi
wL=ws2;   wH=pi;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{3} = [1,0;0,-1];
prob.Psi{3} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{3} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{3} = [1;zeros(n-1,1)];
Mt = sparse([zeros(n,n+1); zeros(1,n) -1]);
prob.M{3} = [Mn Mt(:)];
prob.N{3} = sparse(zeros(n+1,n+1));

% H(w) >= 1-tp,     wp1 <= w <= wp2
wL=wp1;   wH=wp2;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{4} = [1,0;0,-1];
prob.Psi{4} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{4} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{4} = [1;zeros(n-1,1)];
Mt = sparse(zeros(n+1,n+1));
prob.M{4} = [-Mn Mt(:)];
prob.N{4} = sparse([zeros(n,n+1); zeros(1,n), 1-tp]);

% H(w) <= 1+tp,     wp1 <= w <= wp2
wL=wp1;   wH=wp2;
wC=(wL+wH)/2;   wD=(wH-wL)/2;
prob.Phi{5} = [1,0;0,-1];
prob.Psi{5} = [0,exp(j*wC);exp(-j*wC),-2*cos(wD)];
prob.A{5} = [zeros(1,n);eye(n-1),zeros(n-1,1)];
prob.B{5} = [1;zeros(n-1,1)];
Mt = sparse(zeros(n+1,n+1));
prob.M{5} = [Mn Mt(:)];
prob.N{5} = sparse([zeros(n,n+1); zeros(1,n), -1-tp]);

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
ts = sol.x(n+2);

wN = 1000;
wr = linspace(0,pi,wN).';
Hwr = repmat(x0,wN,1) + 2*real((repmat(exp(-j*wr),1,n).^repmat([1:n],wN,1))*x1n);

% plot filter gain
figure;
semilogy([0,ws1,ws1],[ts,ts,ts*0.1],'r--',...
    [ws2,ws2,pi],[ts*0.1,ts,ts],'r--',...
    [wp1,wp1,wp2,wp2,wp1],[1-tp,1+tp,1+tp,1-tp,1-tp],'r--');
hold on;
semilogy(wr,abs(Hwr));
hold off;
xlabel('w (radian)');
ylabel('|H(w)|');
axis([0,pi,ts*0.1,2]);
