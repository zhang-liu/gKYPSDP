function [sos_data,err] = to_sos(prob, sample, feasx)
%
% function [sos_data,err] = to_sos(prob, sample, feasx);
%
% Reformulate discrete-time KYP LMI into SOS sampling version

L=prob.L;
n=prob.n;
m=prob.m;
p=prob.p;
nm=prob.nm;

sos_data.qs=prob.w;
sos_data.varargin.y = feasx;

for i=1:L
    if n{i}~=0
        % compute samples
        
        if sample == 1 | sample == 2
            if max(abs(eig(prob.A{i}))) > 0.9999999
                disp(sprintf('A is not stable by LQR after transform to discrete, max eig = %f',max(abs(eig(prob.A{i})))));
                % disp('switch to uniform sampling on the unit circle');
            end
            A=prob.A{i};
            B=prob.B{i}(:,1);
            R=dlyapchol(A,B);
            Rinv=inv(R);
            A=Rinv'*A*R';
            B=Rinv'*B;
            C = B'*(eye(n{i})+A')^-1*(eye(n{i})+A);
            D = B'*(eye(n{i})+A')^-1*B - eye(1);
            if sum(sum(abs([A,B;C,D]*[A,B;C,D]'-eye(n{i}+1))))>1E-2 | cond(R) > 1E15
                if sum(sum(abs([A,B;C,D]*[A,B;C,D]'-eye(n{i}+1))))>1E-2
                    disp(sprintf('A,B is not balanced, res = %f',sum(sum(abs([A,B;C,D]*[A,B;C,D]'-eye(n{i}+1))))));
                end
                if cond(R) > 1E15
                    disp(sprintf('controllability grammian is not positive definite, cond(R) = %f',cond(R)));
                end
                % disp(sprintf('switch to uniform sampling on the unit circle,i=%d',i));                
            end
            prob.Kt{i}=prob.Kt{i}*[R',zeros(n{i},m{i});zeros(m{i},n{i}),eye(m{i},m{i})];
            if sample == 1 % best method of sampling, balanced SS
                Nz = [A',C'*D,C'*C; B',D'*D,D'*C; zeros(n{i},n{i}),B,A];
                %[X,Z] = eig(Nz);
                [X,Z] = schur(Nz);
                [X,Z] = rsf2csf(X,Z);
                z = diag(Z);
                u = X(n{i}+1,:);
                S = X(n{i}+2:2*n{i}+1,:)*diag(1./u);
            elseif sample == 2 % uniform sampling, balanced SS
                ns = 2*n{i}+1;
                z = linspace(-pi,pi,ns+2).';
                z = exp(j*z(2:end-1,1));
                S = zeros(n{i},ns);
                for ii=1:ns
                    S(:,ii)=(z(ii,1)*eye(n{i})-A)\B;
                end
            end
        elseif sample == 3 % uniform sampling, original SS
            ns = 2*n{i}+1;
            z = linspace(-pi,pi,ns+2).';
            z = exp(j*z(2:end-1,1));
            S = zeros(n{i},ns);
            A = prob.A{i};
            B = prob.B{i}(:,1);
            for ii=1:ns
                S(:,ii)=(z(ii,1)*eye(n{i})-A)\B;
            end    
        end
        
        % transform the first part for kyp
        ns = size(S,2);
        if m{i}==1
            Wu = [S; ones(1,ns)];
        elseif m{i}==2
            Wu = [S; ones(1,ns); zeros(1,ns)];
        end
        
        KtWu = prob.Kt{i}*Wu;
        Wn = diag(KtWu'*prob.N{i}*KtWu);
        Wm = zeros(ns,p);
        for ii =1:p
            Wm(:,ii) = sum(KtWu'*reshape(prob.M{i}(:,ii),nm{i},nm{i}).*(KtWu.'),2);
        end
        WAs = eye(ns);
        
        if m{i} == 2
            V2 = [zeros(n{i}+1,1);1];
            V1 = [eye(n{i}+1);zeros(1,n{i}+1)];
            %V1 = Wu;
            Wm2 = zeros(1,p);
            Wm12 = zeros(size(V1,2),p);
            KtV2 = prob.Kt{i}*V2;
            KtV1 = prob.Kt{i}*V1;
            Wn2 = diag(KtV2'*prob.N{i}*KtV2);
            Wn12 = diag(KtV1'*prob.N{i}*repmat(KtV2,1,size(V1,2)));
            for ii=1:p
                Wm2(:,ii) = diag(KtV2'*reshape(prob.M{i}(:,ii),nm{i},nm{i})*KtV2);
                Wm12(:,ii) = diag(KtV1'*reshape(prob.M{i}(:,ii),nm{i},nm{i})*repmat(KtV2,1,size(V1,2)));
            end
            Wu = [Wu,V2,V1+repmat(V2,1,size(V1,2)),V1-repmat(V2,1,size(V1,2)),...
                V1+repmat(V2,1,size(V1,2))*j,V1-repmat(V2,1,size(V1,2))*j];
            Wn = [Wn;Wn2;Wn12;zeros(size(Wn12))];
            Wm = [Wm;Wm2;Wm12;zeros(size(Wm12))];
            WAs = [eye(ns+1), zeros(ns+1,4*(size(V1,2))); ...
                zeros(size(V1,2),ns+1),0.25*eye(size(V1,2)),-0.25*eye(size(V1,2)),zeros(size(V1,2),2*(size(V1,2))); ...
                zeros(size(V1,2),ns+1+2*size(V1,2)),-0.25*eye(size(V1,2)),0.25*eye(size(V1,2))];
        end 
        
        sos_data.Qs{i}{1} = zeros(nm{i},nm{i});
        sos_data.As{i}{1} = WAs;
        sos_data.Cs{i}{1} = Wu';
        sos_data.Bs{i} = real(Wm);
        sos_data.bs{i} = -(real(Wn));
        
        % determine initial starting point     
        epsI = eye(nm{i})*1E-6;  % purpose: strictly feasible
        KtNKt = -prob.Kt{i}'*(prob.N{i}+reshape(prob.M{i}*sos_data.varargin.y,nm{i},nm{i}))*prob.Kt{i};
        epsImN = epsI + KtNKt;
        if m{i} == 2
            B = [B,zeros(n{i},1)];
        end
        [X1,L1,G1,report1] = dare(A,B,epsImN(1:n{i},1:n{i}),epsImN(n{i}+1:nm{i},n{i}+1:nm{i}),epsImN(1:n{i},n{i}+1:nm{i}),eye(n{i}));
        if report1 >= 0
            X2 = [A'*X1*A-X1, A'*X1*B; B'*X1*A, B'*X1*B]+epsImN;
            [V2,D2]=eig(X2);
            ttt = max(diag(D2))*1E-2;
            D2 = diag((diag(D2)<ttt).*ttt+(diag(D2)>=ttt).*diag(D2));
            X2 = V2*D2*V2';
            X2 = 0.5*(X2+X2');
        end
        
        if report1 >= 0 & ttt > 1E-2
            sos_data.varargin.X{i}{1} = X2; 
        else
            if m{i}==1
                sos_data.varargin.X{i}{1} = eye(n{i}+1);
            elseif m{i}==2
                if KtNKt(end,end) < 0.01
                    KtNKt(end,end) = 0.01;
                end  
                tempTT = KtNKt(1:n{i}+1,end)*(KtNKt(end,end)^(-1))*KtNKt(end,1:n{i}+1);
                sos_data.varargin.X{i}{1} = [eye(n{i}+1)*max([max(eig(tempTT))+0.01,1]),KtNKt(1:n{i}+1,end);KtNKt(end,1:n{i}+1),KtNKt(end,end)];
            end
        end
        
        % transform the second part for gkyp
        if sum(sum(abs(prob.Psi{i}))) > 10E-6
            Wv = S;   
            wtil = real(prob.Psi{i}(1,1)*(conj(z).*z)+prob.Psi{i}(1,2)*conj(z)...
                +prob.Psi{i}(2,1)*z+prob.Psi{i}(2,2)*ones(size(z)));
            WAs = diag(wtil);
            
            if m{i} == 2
                WAs = [WAs;zeros(size(Wn,1)-size(WAs,1),size(WAs,2))];
            end
            
            sos_data.Qs{i}{2} = zeros(n{i},n{i});
            sos_data.As{i}{2} = WAs;
            sos_data.Cs{i}{2} = Wv';
            sos_data.varargin.X{i}{2} = eye(n{i});
        end
    end
end
err=0;

    