function ex_PID()
%
% function ex_PID();
%
% PID controller design
% User can specify different controller design parameters in ex_PID.m
% Refer to Section 4.5 in the gKYPSDP user guide.

disp('**************** PID controller design *********************');

% controller design paramters
wl = 0.4;
wh = 4.0;
th = 0.1;
a = -4;
b = 1;
c = -1;
Ts = 0.1;

th2 = th^2;

% Plant SS
Ps = tf([10],conv(conv([1,1],[1,1,10]),[1,4,15]));
Pss = ss(Ps);
Ap = Pss.A;
Bp = Pss.B;
Cp = Pss.C;
Dp = Pss.D;

% Controller SS
Ak = [0,1;0,-1/Ts];
Bk = [0;1];
% Ck = [ki/Ts,ki-kd/Ts^2];
% Dk = [kp+kd/Ts];
Ck = [0,0; 1/Ts,1; 0,-1/Ts^2];
Dk = [1; 0; 1/Ts];

% Open loop SS
Al = [Ap,zeros(size(Ap,1),size(Ak,2));Bk*Cp,Ak];
Bl = [Bp;Bk*Dp];
%Cl = [Dk*Cp,Ck];
%Dl = [Dk*Dp];
Cl = [Dk*Cp,Ck];
Dl = [Dk*Dp];
n = size(Al,1);

% construct gKYP-SDP
%prob.p = 4;
prob.L = 3;
%prob.n = {n,n,n};
%prob.m = {2,1,1};
prob.w = [0;0;0;-1]; % x = [kp;ki;kd;tl]

% |L(w)| <= th,        wh <= w <= \infty
prob.Phi{1} = [0,1;1,0];
prob.Psi{1} = [0,j;-j,-2*wh];
prob.A{1} = Al;
prob.B{1} = [Bl,zeros(n,1)];
for jj = 1:3
    Mt = sparse([zeros(n,n),zeros(n,1),Cl(jj,:)';zeros(1,n),0,Dl(jj,:)';Cl(jj,:),Dl(jj,:),0]);
    prob.M{1}(:,jj) = Mt(:);
end
Mt = sparse(zeros(n+2,n+2));
prob.M{1}(:,jj+1) = Mt(:);
prob.N{1} = sparse(diag([zeros(n,1);-th2;-1]));

% a Re(L(w)) + b Im(L(w)) + c <= 0,     0 <= w <= \infty
prob.Phi{2} = [0,1;1,0];
prob.Psi{2} = [0,j;-j,0];
prob.A{2} = Al;
prob.B{2} = Bl;
at = a; bt = b; ct = c;
for jj = 1:3
    Mt = sparse([zeros(n,n),Cl(jj,:)'*(at+bt*j); (at-bt*j)*Cl(jj,:),(at+j*bt)*Dl(jj,:)'+(at-j*bt)*Dl(jj,:)]);
    prob.M{2}(:,jj) = Mt(:);
end
Mt = sparse(zeros(n+1,n+1));
prob.M{2}(:,jj+1) = Mt(:);
prob.N{2} = sparse([zeros(n,n+1);zeros(1,n),2*ct]);

% Im(L(w)) <= -tl,      0 <= w <= wl
prob.Phi{3} = [0,1;1,0];
prob.Psi{3} = [-1,j*wl/2;-j*wl/2,0];
prob.A{3} = Al;
prob.B{3} = Bl;
at = 0; bt = 1;
for jj = 1:3
    Mt = sparse([zeros(n,n),Cl(jj,:)'*(at+bt*j); (at-bt*j)*Cl(jj,:),(at+j*bt)*Dl(jj,:)'+(at-j*bt)*Dl(jj,:)]);
    prob.M{3}(:,jj) = Mt(:);
end
Mt = sparse([zeros(n,n+1);zeros(1,n),2]);
prob.M{3}(:,jj+1)=Mt(:);
prob.N{3} = sparse(zeros(n+1,n+1));

if 1==1
    % special-purpose solver, gkypsdp
    [sol,info] = gkypsdp_solver(prob);
else
    % general-purpose solver, sedumi
    opt.IPMsolver = 2;
    [sol,info] = gkypsdp_solver(prob,opt);
end

kp = sol.x(1);
ki = sol.x(2);
kd = sol.x(3);
tl = sol.x(4);

Cl = [kp,ki,kd]*Cl;
Dl = [kp,ki,kd]*Dl;
SSl = ss(Al,Bl,Cl,Dl);
SScl = feedback(SSl,1);

wN = 500;
ws = logspace(log10(wl/2),log10(10*wh),wN)';
Lws = zeros(wN,1);
for ii = 1:wN
    Lws(ii,1) = Cl*((j*ws(ii,1)*eye(n)-Al)\Bl)+Dl;
end

% Nyquist plot
figure;
plot(Lws);
hold on;
plot([Cl*((j*wl*eye(n)-Al)\Bl)+Dl;Cl*((j*wh*eye(n)-Al)\Bl)+Dl],'o');
plot(th*exp(j*linspace(0,2*pi,50)),'r');
plot([-2+((-a/b)*(-2)-c/b)*j; 1+((-a/b)*(1)-c/b)*j],'r');
plot([-2+(-tl)*j; 1+(-tl)*j],'r');
plot([-1+0.000001*j],'*');plot([-1+0.000001*j],'o');
axis([-1.5,0.5,-1.5,0.5]);
grid on;
hold off;
xlabel('Re(L(w))');
ylabel('Im(L(w))');

% Plot step response
figure;
step(SScl);

