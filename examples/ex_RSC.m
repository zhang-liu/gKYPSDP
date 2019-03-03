function ex_RSC()
%
% function ex_RSC();
%
% Robust stabilizing controller (RSC) synthesis
% Simultaneous stabilization of 2^4 = 16 plant vertices of ARMAX model of
% PUMA 762 robotic disk grinding process.
% Refer to Section 4.6 in the gKYPSDP user guide.

disp('********* Robust stabilizing controller (RSC) synthesis *********');

% Vertex systems
nv = 2^4;
S = cell(nv,1);
k = 1;
q10 = 0.0257; q20 = -0.0764; q30 = -0.1619; q40 = -0.1688;
error = 20; % in percent around nominal values
for q1 = q10+abs(q10)*[-1 1]*error/100
  for q2 =  q20+abs(q20)*[-1 1]*error/100
    for q3 = q30+abs(q30)*[-1 1]*error/100
      for q4 = q40+abs(q40)*[-1 1]*error/100
	S{k} = tf([q1 q2 q3 q4 0 0 0 0],conv([1 -1],...
                  [1 -1.914 1.779 -1.0265 0.2508 0 0 0 0 0 0 0]),1);
        k = k+1;
      end
    end
  end
end

% Central polynomial
c = [zeros(1,19) 1];

% Stability region
H = [1 0;0 -1];

% Construct gKYP-SDP
strict = 1e-3;
prob.L = nv;
%prob.p = 7+8; % [x,y]
n = 19;
%prob.n = {n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n};
%prob.m = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
prob.w = zeros(7+8,1); 
AB = [zeros(19,1) eye(19)];
II = eye(16);
for k = 1:nv
    prob.A{k} = AB(1:19,1:19);
    prob.B{k} = AB(1:19,end);
    prob.Phi{k} = H;
    prob.Psi{k} = zeros(2,2);
    a = get(S{k},'den'); a = fliplr(a{1});
    b = get(S{k},'num'); b = fliplr(b{1});
    % Build Sylvester matrix with monic controller denominator
    T = [toeplitz([a(1);zeros(7,1)],[a zeros(1,7)]);
         toeplitz([b(1);zeros(7,1)],[b zeros(1,7)])];
    for l = 1:7
        Mt = sparse(-c'*II(l,:)*T-(c'*II(l,:)*T)');
        prob.M{k}(:,l) = Mt(:);
    end
    l = 8;
    prob.N{k} = sparse(strict*eye(20)-c'*II(l,:)*T-(c'*II(l,:)*T)');
    for l = 9:16
        Mt = sparse(-c'*II(l,:)*T-(c'*II(l,:)*T)');
        prob.M{k}(:,l-1) = Mt(:);
    end
end

if 1==1
    % special-purpose solver, gkypsdp
    [sol,info] = gkypsdp_solver(prob);
else
    % general-purpose solver, sedumi
    opt.IPMsolver = 2;
    [sol,info] = gkypsdp_solver(prob,opt);
end

x = sol.x(1:7)';
y = sol.x(8:end)';
C = tf(fliplr(y),fliplr([x 1]),1);

% Plot random robust root locii in black
figure;
th = 0:0.01:2*pi;
plot(cos(th),sin(th),'k--');     
axis equal; hold on;
for k = 1:5,
    for i = 1:100,
        q1 = q10+abs(q10)*(-1+2*rand)*error/100;
        q2 = q20+abs(q20)*(-1+2*rand)*error/100;
        q3 = q30+abs(q30)*(-1+2*rand)*error/100;
        q4 = q40+abs(q40)*(-1+2*rand)*error/100;
        SC = feedback(tf([q1 q2 q3 q4 0 0 0 0], conv([1 -1],...
            [1 -1.914 1.779 -1.0265 0.2508 0 0 0 0 0 0 0]), 1), C);
        [z,p] = zpkdata(SC); p = p{1};
        h = plot(real(p),imag(p),'.k'); set(h,'MarkerSize',10);
    end
    drawnow
end
% Vertex root locii in red
for k = 1:nv
    SC = feedback(S{k},C);
    [z,p,k] = zpkdata(SC); p = p{1};
    h = plot(real(p),imag(p),'.r'); set(h,'MarkerSize',10);
end
% Nominal root locii in blue
SC = feedback(tf([q10 q20 q30 q40 0 0 0 0], conv([1 -1],...
            [1 -1.914 1.779 -1.0265 0.2508 0 0 0 0 0 0 0]), 1), C);
[z,p] = zpkdata(SC); p = p{1};
h = plot(real(p),imag(p),'.b'); set(h,'MarkerSize',15);
xlabel('Real');  ylabel('Imaginary');
hold off;
