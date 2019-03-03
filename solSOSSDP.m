function [y,X,z,ptime,stime,iters,p_feas,d_feas,abs_gap,rel_gap] = solSOSSDP(A,B,b,C,q,Q,options,varargin)
%
% [y,X,z,ptime,stime,iters,p_feas,d_feas,abs_gap,rel_gap] =
% solSOSSDP(A,B,b,C,q,Q,options,varargin);
%
% solves the SDP
%
% Primal
%    minimize    q'*y + sum_i sum_k tr(Q{i}{k}*X{i}{k})
%    subject to  sum_k(i) A{i}{k}*diag(C{i}{k}*X{i}{k}*C{i}{k}') 
%                    + B{i}*y = b{i},    i=1,...,L
%                X{i}{k} >= 0,   i=1,...,L,  k=1,...,k(i)
%
% Dual
%    maximize    sum_i b{i}'*z{i}
%    subject to  C{i}{k}'*diag(A{i}{k}'*z{i})*C{i}{k} <= Q{i}{k},
%                    i=1,...,L,  k=1,...,k(i)
%                sum_i B{i}'*z{i} = q
%
% using infeasible primal-dual interior-point method.

%
% Input:
%    constr.A{i}{k}: m(i)-by-l(i), real
%    constr.B{i}: m(i)-by-p, real
%    constr.b{i}: m(i)-by-1, real
%    constr.C{i}{k}: l(i)-by-n{i}(k), complex
%    q: p-by-1, real
%    Q{i}{k}: n{i}(k)-by-n{i}(k), complex
%
% Output:
%    y: primal (linear) solution, real
%    X{i}{k}: primal (SDP) solution, hermitian
%    z{i}: dual solution, real
%    ptime: preprocessing time
%    stime: solving time
%    iters: total number of iterations

%warning off

%%% PARAMETERS %%%
ABSTOL = options.abstol; %1e-6;
RELTOL = options.reltol; %1e-6;
PINFEASTOL = options.pfeastol; %1e-6;
DINFEASTOL = options.dfeastol; %1e-6;
RCONDTOL = 1E-18;
RCONDTOLITERS = 12; 
MAXITERS = options.maxiters; %50
STEP = 0.9;
EXPON = 3;

NUMTEST = 3; % set this number to > MAXITERS, then disable feasibility check
GAPLIMIT = 4^NUMTEST;
obj_gap = ones(NUMTEST,1);
prim_div = zeros(NUMTEST,1);
dual_div = zeros(NUMTEST,1);

m = cellfun('size',B,1)';
p = length(q);
L = length(m);
for i=1:L
   K(i) = length(C{i});
   l{i} = cellfun('size',C{i},1);
   n{i} = cellfun('size',C{i},2);
end

%%% PREPROCESSING STEPS %%%
ptime = cputime;

prob_size = 0;
for i=1:L
   prob_size = prob_size + sum(n{i});
end

% Pre-compute norm of b, Q
for i=1:L
   norm_b(i) = norm(b{i});
   for k=1:K(i)
      norm_Q(i,k) = norm(full(Q{i}{k}));
   end
end
norm_q = norm(q);


ptime = cputime - ptime;

%%% DETERMINE INITIAL POINTS %%%
if (nargin > 7) % add code to check that x,P,Z are feasible
   y = varargin{1}.y;  
   X = varargin{1}.X;
   if isfield(varargin,'z') == 1
      z = varargin.z;
   else
      for i=1:L
         z{i} = zeros(m(i),1);
      end
   end
   
   for i=1:L
      for k=1:K(i)
         %S{i}{k} = Q{i}{k} - C{i}{k}'*diag(A{i}{k}'*z{i})*C{i}{k};
         S{i}{k} = eye(n{i}(k));
      end
   end
else 
   y = zeros(p,1);
   for i=1:L
      z{i} = zeros(m(i),1);
      for k=1:K(i)
         X{i}{k} = eye(n{i}(k));
         S{i}{k} = eye(n{i}(k));
      end
   end
end

disp([' ']);
disp(['#ITERS   PRIMAL INF.    DUAL INF.      PRIMAL OBJ.       DUAL OBJ.      REL. GAP       ABS. GAP    CPU TIME']); 
disp([' ']);


startime = cputime;
mu=NaN;
prev_relprimalres = Inf;
prev_dualres  = Inf;

for iters=0:MAXITERS

   % STOPPING CRITERIA
   primalcost = 0;
   dualcost = -q'*y;
   gap = 0;
   primalres2 = q;
   for i=1:L
      primalcost = primalcost - b{i}'*z{i};
      primalres2 = primalres2 - B{i}'*z{i};
      dualres{i} = b{i} - B{i}*y;
      for k=1:K(i)
         dualcost = dualcost - Q{i}{k}(:)'*X{i}{k}(:);
         gap = gap + S{i}{k}(:)'*X{i}{k}(:);
         primalres1{i}{k} = ...
            Q{i}{k} - C{i}{k}'*(diag(A{i}{k}'*z{i}))*C{i}{k} - S{i}{k};
         norm_primalres1(i,k) = norm(primalres1{i}{k});
         dualres{i} = dualres{i} - ...
            A{i}{k}*diag(C{i}{k}*X{i}{k}*C{i}{k}');
      end
      norm_dualres(i) = norm(dualres{i});
   end
   abs_gap(iters+1) = gap;
   
   if (primalcost < 0)  
      relgap = gap/(-primalcost);
   elseif (dualcost > 0)
      relgap = gap/(dualcost);
   else
      relgap = Inf;
   end
   rel_gap(iters+1) = relgap;

   relprimalres = norm(primalres2)/max(1,norm_q);
   reldualres = 0;
   for i=1:L, 
      for k=1:K(i)
         relprimalres = max( relprimalres, ...
            norm(primalres1{i}{k})/max(1,norm_Q(i,k)) );
      end;
      reldualres = max(reldualres, norm(dualres{i})/max(1,norm_b(i)));
   end;

%   if ((relprimalres > prev_relprimalres) | (reldualres > prev_dualres))
%       divflag = 1;  
%   else, divflag = 0;
%   end;
   prev_relprimalres = relprimalres;
   prev_reldualres = reldualres;
   
   p_feas(iters+1) = relprimalres;
   d_feas(iters+1) = reldualres;
  
   stime = cputime - startime;
   
      disp([sprintf('% 4.0d',  iters), ...
         sprintf('% 16.5e', relprimalres), ...
         sprintf('% 15.5e', reldualres), ...
         sprintf('% 16.5e', primalcost), ...
         sprintf('% 16.5e', dualcost), ...
         sprintf('% 15.5e', relgap), ...
         sprintf('% 15.5e', gap), ...
         sprintf('% 8.1f',  stime)]);

   %%%%%%%%%%%%%%%%%%%%% infeasible detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   obj_gap(1:end-1,1) = obj_gap(2:end,1);
   obj_gap(NUMTEST,1) = abs(primalcost-dualcost);
   prim_div(1:end-1,1) = prim_div(2:end,1);
   prim_div(NUMTEST,1) = relprimalres;
   dual_div(1:end-1,1) = dual_div(2:end,1);
   dual_div(NUMTEST,1) = reldualres;
   if iters > RCONDTOLITERS & sum(diff(obj_gap)>0) == NUMTEST-1  & obj_gap(end,1)/obj_gap(1,1)>GAPLIMIT
       disp('Numerical problem: duality gap is diverging.');  break;
   end
   if iters > RCONDTOLITERS & sum(diff(prim_div)>0) == NUMTEST-1 & prim_div(end,1) > min([gap,relgap]) % prim_div(end,1) > 1E-4
       disp('Numerical problem: primal infeasibility residue is diverging'); break;
   end
   if iters > RCONDTOLITERS & sum(diff(dual_div)>0) == NUMTEST-1 & dual_div(end,1) > min([gap,relgap]) % dual_div(end,1) > 1E-4
       disp('Numerical problem: dual infeasibility residue is diverging'); break;
   end
    
   if ((relprimalres <= PINFEASTOL) & (reldualres <= DINFEASTOL))
      if (gap <= ABSTOL)
         disp('Absolute tolerance reached.');  break;
      elseif (relgap <= RELTOL)
         disp('Relative tolerance reached.');  break;
      end;
   end;
   %if dualcost > 0, break; end  %% ??? %%
   %if (relgap <= max(relprimalres, reldualres)),
   %   disp('Relative gap less than relative residuals.'); 
   %   stop=0;
   %   keyboard;  if (stop), return; end;
   %end;

   Bsc = zeros(sum(m),p);  d3sc = zeros(sum(m),1);
   indBsc=0;

   for i=1:L

      d3{i} = dualres{i};
      Atsc = zeros(sum(l{i}),m(i));
      indAtsc = 0;

      % SCALING: compute R, lambda
      for k=1:K(i)

         [L1,flag] = chol(S{i}{k});
         if (flag), 
            disp(['Numerical problem: S{' num2str(i) ',' num2str(k) '} is not positive definite.']); 
            %keyboard; 
            return; 
         end;
         L1 = L1';  % Sik = L1*L1'

         [L2,flag] = chol(X{i}{k});
         if (flag), 
            disp(['Numerical problem: X{' num2str(i) ',' num2str(k) '} is not positive definite.']); 
            %keyboard; 
            return; 
         end;
         L2 = L2';  % Xik = L2*L2'

         % R'*Z*R = Rinv*S*Rinv' = lambda
         [UU,Sig,VV] = svd(L2'*L1);
         lambda{i}{k} = diag(Sig);
         if (min(lambda{i}{k}) < 0), 
            disp('Lambda negative.'); 
            keyboard;
         end;
         R{i}{k} = L1*VV*diag(1./sqrt(lambda{i}{k}));
         Rinv{i}{k} = diag(1./sqrt(lambda{i}{k}))*UU'*L2';

         % scaled data
         % - Csc = C*Rinv'
         % - Qsc = Rinv*Q*Rinv'
         Csc{i}{k} = C{i}{k}*Rinv{i}{k}'; 
         Qsc{i}{k} = Rinv{i}{k}*Q{i}{k}*Rinv{i}{k}';  
                % debug for problems with nonzero Q
         primalres1sc{i}{k} = Rinv{i}{k}*primalres1{i}{k}*Rinv{i}{k}';
              
         % FACTOR COEFFICIENT MATRIX
         %
         % We will need to solve two equations of the form
         %
         %    A*diag(C*T*C'*diag(A'*dz)*C*T*C') + B*dy = d3
         %                                       B'*dz = d2
         %   
         % solve [H{1}   0     0   ...    0    B{1}][dz{1}]   [d3{1}]
         %       [ 0    H{2}   0   ...    0    B{2}][dz{2}]   [d3{2}]
         %       [ 0     0    H{3} ...    0    B{3}][dz{3}] = [d3{3}]
         %       [ :     :     :          :     :  ][ :   ]   [  :  ]
         %       [ 0     0     0         H{L}  B{L}][dz{L}]   [d3{L}]
         %       [B{1}' B{2}' B{3}' ...  B{L}'   0 ][ dy  ]   [ d2  ]
         %
         
         %%%%%%%%%%%%%%%%%%%% sqr for complex hermitian matrix %%%%%%%%%%%%%%
         % sqrCsc = (Csc{i}{k}*Csc{i}{k}').^2;
         sqrCsc = abs(Csc{i}{k}*Csc{i}{k}').^2;
         if max(max(isnan(sqrCsc))) > 0
             break;
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % CR*CR.^2 = VCR*VCR'
         [VCR,DCR] = eig(sqrCsc);
         VCR = VCR*diag(sqrt(max(0,diag(DCR))));      
         Atsc(indAtsc+[1:l{i}(k)],:) = VCR'*A{i}{k}';
         indAtsc=indAtsc+l{i}(k);


         % AFFINE SCALING EQUATIONS
         % 
         %    -Tinv*dXa*Tinv + C'*diag(A'*dza)*C = rhs1
         %              A*diag(C*dXa*C') + B*dya = rhs2
         %                                B'*dza = rhs3
         %
         % where rhs1 = primalres1 - 2*R*(-diag(lambda)^2 .* G)*R' 
         %              (where G_ij = 1/(lambdai + lambdaj) )
         %            = primalres1 + R*diag(lambda)*R'
         %            = primalres1 + S
         %            = Q - C'*diag(A'*z)*C 
         %       rhs2 = dualres
         %       rhs3 = primalres2
         % 
         % Solve by eliminating dXa:
         %
         %   [ H  B ][dza] = [ rhs2 + A*diag(C*T*rhs1*T*C') ] = [ d3 ]
         %   [ B' 0 ][dya]   [         rhs3                 ]   [rhs3] 
         %
         % eliminate dza
         %
         %   B'*H^{-1}*B * dya = B'*H^{-1}*d3 - rhs3
         % 
         %  dza = H^{-1}*(d3 - B*dya)
         %  dXa = T*C'*diag*(A'*dza)*C*T - T*rhs1*T
         

         % rhs1 = Rinv * (primalres1 + S) * Rinv'
         rhs1{i}{k} = Qsc{i}{k} - ...
                 Csc{i}{k}'*diag(A{i}{k}'*z{i})*Csc{i}{k};  
         d3{i} = d3{i} + A{i}{k}*diag(Csc{i}{k}*rhs1{i}{k}*Csc{i}{k}');
         
      end

      % B'*H^{-1}*B = LH'*LH
      LH{i} = triu(qr(Atsc,0));  LH{i} = LH{i}(1:m(i),1:m(i));
      if rcond(LH{i}) < RCONDTOL & iters > RCONDTOLITERS
          disp('Numerical problem: rcondLH.');
          return;
      end
      Bsc(indBsc+[1:m(i)],:) = LH{i}'\B{i}; 
      d3sc(indBsc+[1:m(i)]) = LH{i}'\d3{i};
      indBsc=indBsc+m(i);

   end

   % B'*H^{-1}*B = LHH*LHH'
   [QHH,LHH] = qr(Bsc,0);  QHH= QHH(:,1:p); LHH=LHH(1:p,1:p);
   
   dya = -LHH\(LHH'\primalres2 - QHH'*d3sc);
   indBsc=0;

   for i=1:L

      dza{i} = LH{i} \ (d3sc(indBsc+[1:m(i)]) - ...
                        Bsc(indBsc+[1:m(i)],:)*dya);
      %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
      dza{i} = real(dza{i});
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      indBsc=indBsc+m(i);

      for k=1:K(i)

         CAtzC = C{i}{k}'*diag(A{i}{k}'*dza{i})*C{i}{k};
         CscAtzCsc = Csc{i}{k}'*diag(A{i}{k}'*dza{i})*Csc{i}{k};
         %CscAtzCsc = Rinv{i}{k}*CAtzC*Rinv{i}{k}';

         dSa{i}{k} = primalres1{i}{k} - CAtzC;
         dSa{i}{k} = 0.5*(dSa{i}{k} + dSa{i}{k}');
         dS{i}{k} = primalres1sc{i}{k} - CscAtzCsc;
         dS{i}{k} = 0.5*(dS{i}{k} + dS{i}{k}');

         dXa{i}{k} = -Rinv{i}{k}'*(rhs1{i}{k} - CscAtzCsc)*Rinv{i}{k};
         dXa{i}{k} = 0.5*(dXa{i}{k} + dXa{i}{k}');
         dX{i}{k} = -(rhs1{i}{k} - CscAtzCsc);
         dX{i}{k} = 0.5*(dX{i}{k} + dX{i}{k}');

         %eig_X(i,k) = max(-real(eig(dXa{i}{k},X{i}{k})));
         eig_X(i,k) = max(-real(eig(dX{i}{k},diag(lambda{i}{k}))));
         %eig_S(i,k) = max(-real(eig(dSa{i}{k},S{i}{k})));
         eig_S(i,k) = max(-real(eig(dS{i}{k},diag(lambda{i}{k}))));

      end

   end

   % STEP LENGTHS TO BOUNDARY with  TTT heuristic
   maxeigp = max(0, max(max(eig_S)));
   maxeigd = max(0, max(max(eig_X)));
   prstep = 1/max(1.0, maxeigp);
   dustep = 1/max(1.0, maxeigd);
   newgap = 0;
   newgap2 = 0;
   for i=1:L
      for k=1:K(i)
%         newgap = newgap + ...
%            (S{i}{k}(:) + prstep*dSa{i}{k}(:))' * ...
%            (X{i}{k}(:) + dustep*dXa{i}{k}(:));
         newgap = newgap + ...
            lambda{i}{k}'*lambda{i}{k} +  ...
            prstep*diag(dS{i}{k})'*lambda{i}{k} + ...
            dustep*lambda{i}{k}'*diag(dX{i}{k}) + ...
            prstep*dustep*dS{i}{k}(:)'*dX{i}{k}(:);
      end
   end

%   if (gap/prob_size > 1e-6),
%      minstep = min(dustep,prstep);
%      delta = max(EXPON, 3*minstep^2)
%   else
      delta = EXPON;
%   end;
 
   mu = (gap/prob_size)*min(1,(newgap/gap)^delta);


   % CENTERING-CORRECTOR EQUATIONS 
   %
   %    -Tinv*dZ*Tinv' + C'*diag(A'*dz)*C = rhs1
   %               A*diag(C*dX*C') + B*dy = rhs2
   %                                B'*dz = rhs3
   %
   % where rhs1 = 0 - 2*R*((mu*I-Hc).*G)*R' 
   %     and Hc = 0.5*(R'*(dSa*dZa)*Rinv' + Rinv*(dZa*dSa)*R) 
   %                  G_ij = 1/(lambda_i + lambda_j)
   %       rhs2 = 0
   %       rhs3 = 0
   %
   % Solve by eliminating dX:
   %
   %   [ H  B ] [dz] = [ A*diag(C*T*rhs1*T*C') ] = [ d3  ]
   %   [ B' 0 ] [dy]   [          0            ]   [  0  ] 
   %
   % eliminate dz
   %
   %   B'*H^{-1}*B * dy = B'*H^{-1}*d3 - rhs3
   % 
   %  dz = H^{-1}*(d3 - B*dy)
   %  dX = T*C'*diag*(A'*dz)*C*T - T*rhs1*T
   %
   
   indd3sc = 0;
   for i=1:L
      d3{i} = zeros(m(i),1);
      for k=1:K(i)
%         Hc = 0.5*(R{i}{k}'*dXa{i}{k}*dSa{i}{k}*Rinv{i}{k}' + ...
%                   Rinv{i}{k}*dSa{i}{k}*dXa{i}{k}*R{i}{k});
         Hc = 0.5*(dX{i}{k}*dS{i}{k} + dS{i}{k}*dX{i}{k});
         rhs1{i}{k} = -2*((mu*eye(n{i}(k)) - Hc) ./ ...
                  (lambda{i}{k}(:,ones(1,n{i}(k))) + ...
                   lambda{i}{k}(:,ones(1,n{i}(k)))'));  
         d3{i} = d3{i} + A{i}{k}*diag(Csc{i}{k}*rhs1{i}{k}*Csc{i}{k}');
      end
      d3sc(indd3sc+[1:m(i)]) =  LH{i}' \ d3{i};
      indd3sc = indd3sc+m(i);
   end

   dy= LHH\(QHH'*d3sc);
   %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
   dy = real(dy);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   indBsc=0;
   for i=1:L

      dz{i} = LH{i}\(d3sc(indBsc+[1:m(i)]) - Bsc(indBsc+[1:m(i)],:)*dy);
      %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
      dz{i} = real(dz{i});
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      indBsc = indBsc+m(i);

      for k=1:K(i)

         CAtzC = C{i}{k}'*diag(A{i}{k}'*dz{i})*C{i}{k};
         CAtzC = .5*(CAtzC + CAtzC');
         CscAtzCsc = Rinv{i}{k} * CAtzC * Rinv{i}{k}';
         %CscAtzCsc = .5*(CscAtzCsc + CscAtzCsc');
         %CscAtzCsc = Csc{i}{k}'*diag(A{i}{k}'*dz{i})*Csc{i}{k};

         dX{i}{k} = -Rinv{i}{k}'*(rhs1{i}{k} - CscAtzCsc)*Rinv{i}{k};
         dX{i}{k} = dXa{i}{k} + dX{i}{k};
         dX{i}{k} = 0.5*(dX{i}{k} + dX{i}{k}');

         %dX{i}{k} = Rinv{i}{k}' * (dX{i}{k} - rhs1{i}{k} ...
         %          + CscAtzCsc) * Rinv{i}{k};
         %dX{i}{k} = 0.5*(dX{i}{k} + dX{i}{k}');
         eig_X(i,k) = max(-real(eig(dX{i}{k},X{i}{k})));

         dS{i}{k} = -CAtzC;
         dS{i}{k} = 0.5*(dS{i}{k} + dS{i}{k}');
         dS{i}{k} = dSa{i}{k} + dS{i}{k};

         %dS{i}{k} = -CAtzC + R{i}{k}*dS{i}{k}*R{i}{k}';
         %dS{i}{k} = 0.5*(dS{i}{k} + dS{i}{k}');
         eig_S(i,k) = max(-real(eig(dS{i}{k},S{i}{k})));
      end

      dz{i} = dza{i} + dz{i};
      %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
      dz{i} = real(dz{i});
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
   dy = dya + dy;
   %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
   dy = real(dy);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % STEP LENGTHS AND UPDATE
   maxeigp = 0;  maxeigd = 0;
   maxeigp = max(0, max(max(eig_S)));
   maxeigd = max(0, max(max(eig_X)));
   prstep = 1/max(1.0, maxeigp/STEP);
   dustep = 1/max(1.0, maxeigd/STEP);

   y = y + dustep*dy;
   %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
   y = real(y);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   for i=1:L
      z{i} = z{i} + prstep*dz{i};
      %%%%%%%%%%%%%%%%%%%%%  elimiate imaginary part %%%%%%%%%%%%%%%%%%%%%%
      z{i} = real(z{i});
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      for k=1:K(i)
         S{i}{k} = S{i}{k} + prstep*dS{i}{k};
         S{i}{k} = 0.5*(S{i}{k} + S{i}{k}');
         X{i}{k} = X{i}{k} + dustep*dX{i}{k};
         X{i}{k} = 0.5*(X{i}{k} + X{i}{k}');
      end
   end
end

