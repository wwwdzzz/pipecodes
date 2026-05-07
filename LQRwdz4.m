function traj = LQRwdz4(tspan,cont,P,Atil,Btil,datastartmv,datalengthmv,autopart,targetdata,windownum,epsilon,DS,winstart,auto_traj,Q0)

%need to add Qf ,DS.btilde, DS.atilde
%Q = 1*cont.Q(1:size(targetdata,1),1:size(targetdata,1));%%%%%%%0607
n=size(Atil{1},1);
%Q=0.1/50^2*epsilon^2*blkdiag(1*eye(n),0.01*eye(n));
%Q(2,2)=100;
% Q6=0.01*[100 0 0 0 0 0;
%          0 0 0 0 0 0;
%          0 0 100 0 0 0;
%          0 0 0  0 0 0;
%          0 0 0  0 100 0;
%          0 0 0  0 0 0;];
% Q=Q6;
Q=Q0*eye(n);
Rhat =1*(epsilon^2);
%Rhat =0.01*cont.Rhat;
Qf =10*Q;
ep   =epsilon;
%n    = size(DS.M,1);
%tspan =
fullstep=datastartmv{end}+datalengthmv{end}-1;

xautofull=autopart;
xisfull=zeros(size(targetdata,1),20000);

xtfull=autopart;
utfull=zeros(size(Btil{1},2),fullstep);

z0=(targetdata(:,1)-xautofull(:,1))/ep;
%z0=z0(1:size(Atil{1},1),:);

outdof = linspace(1,20,20);



%% The process to get P matrix
%P = [];%%%%%%%%Util
% z0 = Pξ0
for i=winstart:windownum


    %xauto=xautofull(:,datastartmv{i}:fullstep);
    if i==1
        z0temp=z0;%xauto(:,1);
        z0ode=xautofull(:,1);
    else
        z0temp=xtfull(:,datastartmv{i}-1)-xtfull(:,datastartmv{i}-1);%xtfull(:,datastartmv{i}-1)-xautofull(:,datastartmv{i}-1);
        z0ode=xtfull(:,datastartmv{i}-1);
    end
    %xi0   = P{i}'*z0temp;%%%%%%%0607
    xi0   = pinv(P{i})*z0temp;
    Nstep=datalengthmv{i};
    
    utnew=zeros(2,fullstep-datastartmv{i}+1);
    closedtime=linspace(0,size(utnew,2)/100,size(utnew,2))
    prestep=size(utnew,2);
    tpref=size(utnew,2)/100;


    % om = 2*pi/tpref;
    % ufun = @(t) transpose(interp1(closedtime,utnew',t,'linear','extrap'));
    % set(DS,'u',ufun);
    % [tfullauto, xodefullauto] = time_integration_transient(DS,om,'nCycles',...
    %     1, 'nSteps', prestep,...
    %     'integrationMethod','Newmark','outdof',outdof,'init',z0ode);


    %xauto=xautofull;%(1:size(Atil{1},1),datastartmv{i}:fullstep);%%%%%%%%%%%0607
    xauto=xautofull(:,datastartmv{i}:fullstep);
    %xauto=xodefullauto(1:fullstep-datastartmv{i}+1,:)';
   
    %% calculate offline matrix
    bQ_N    = 2*ep*transpose(Qf*P{i})*auto_traj(:,end);%7
    bQ_k    = 2*ep*transpose(Q*P{i})*auto_traj;
    Q2      = ep^2*transpose(P{i})*Q*P{i}; %50x50
    Qfhat   = ep^2*transpose(P{i})*Qf*P{i}; % 50x50

    %Atil = DS.Atil;
    %Btil = DS.Btil;

    M_Nm1    = (Rhat + Btil{i}' * Qfhat * Btil{i})^-1;
    h_Nm1    = -bQ_N' * Btil{i} * M_Nm1 * Btil{i}' * Qfhat' * Atil{i};%0
    K_Nm1    = -Atil{i}' * Qfhat * Btil{i} * M_Nm1 * Btil{i}' * Qfhat' * Atil{i};%1-5
    Qtil_Nm1 = Q2 + Atil{i}' * Qfhat * Atil{i};
    bN_Nm1   = bQ_k(:,end-1)' + bQ_N' * Atil{i};%0


    Mks    = cell(1,fullstep-datastartmv{i}+1);
    hks    = cell(1,fullstep-datastartmv{i}+1);
    Kks    = cell(1,fullstep-datastartmv{i}+1);
    Qtilks = cell(1,fullstep-datastartmv{i}+1);
    bNks   = cell(1,fullstep-datastartmv{i}+1);

    Mks{fullstep-datastartmv{i}+1}    = M_Nm1;
    hks{fullstep-datastartmv{i}+1}    = h_Nm1;
    Kks{fullstep-datastartmv{i}+1}    = K_Nm1;
    Qtilks{fullstep-datastartmv{i}+1} = Qtil_Nm1;
    bNks{fullstep-datastartmv{i}+1}   = bN_Nm1;

    %matlab lqr
    [K,P1,~]=dlqr(Atil{i},Btil{i},Q2,Rhat);



    for k = fullstep-datastartmv{i}:-1:2
        % Mk h(k) Kk Qtilk bNk
        % whether they will converge ?
        % Mks{k}   = (Rhat + Btil{i}' * (Qtilks{k+1} + Kks{k+1}) * Btil{i})^-1;%sl
        % hks{k}   = -(bNks{k+1}+hks{k+1})*Btil{i}*Mks{k}*Btil{i}'*(Qtilks{k+1}+Kks{k+1})'*Atil{i};
        % Kks{k}   = -Atil{i}'*(Qtilks{k+1}+Kks{k+1})*Btil{i}*Mks{k}*Btil{i}'*(Qtilks{k+1}+Kks{k+1})'*Atil{i};
        % Qtilks{k}= Q2 + Atil{i}'*(Qtilks{k+1}+Kks{k+1})*Atil{i};
        % bNks{k}  = bQ_k(:,k)' + (bNks{k+1}+hks{k+1}) * Atil{i};

        Mks{k}   = (Rhat + Btil{i}' * (P1) * Btil{i})^-1;%sl
        hks{k}   = -(bNks{k+1}+hks{k+1})*Btil{i}*Mks{k}*Btil{i}'*(P1)'*Atil{i};
        Kks{k}   = -Atil{i}'*(P1)*Btil{i}*Mks{k}*Btil{i}'*(P1)'*Atil{i};
        Qtilks{k}= Q2 + Atil{i}'*(P1)*Atil{i};
        bNks{k}  = bQ_k(:,k)' + (bNks{k+1}+hks{k+1}) * Atil{i};
    end
    k = k - 1;
    % Mks{k} = (Rhat + Btil{i}' * (Qtilks{2} + Kks{2}) * Btil{i})^-1;   % M0
    Mks{k} = (Rhat + Btil{i}' * (P1) * Btil{i})^-1;
    %% calculate control policy
    us       = zeros(size(Btil{i},2),Nstep);
    xis      = zeros(size(xi0,1), Nstep);
    xis(:,1) = xi0;

    %for j = 1:datalengthmv{i}-1
    for j=1:fullstep-datastartmv{i}
        % us(:,j) = - Mks{j} * ((0.5*(bNks{j+1}+hks{j+1})*Btil{i} + xis(:,j)'*Atil{i}'*(Qtilks{j+1}+Kks{j+1})*Btil{i})');
        us(:,j) = - Mks{j} * ((0.5*(bNks{j+1}+hks{j+1})*Btil{i} + xis(:,j)'*Atil{i}'*(P1)*Btil{i})');
        xis(:,j+1) = Atil{i} * xis(:,j) + Btil{i} * us(:,j);
    end

    j = j + 1;
    us(:,j) = -Mks{j}*((0.5*bQ_N'*Btil{i} + xis(:,j)'*Atil{i}'*Qfhat*Btil{i})');%30



    %% evalute state of original system
    %
    xt     = xauto+ep*real(P{i}*xis);
    xi=real(P{i}*xis);


    for k=1:size(xt,2)
        xtfull(:,datastartmv{i}+k-1)=xt(:,k);
        xisfull(:,datastartmv{i}+k-1)=xi(:,k);
        utfull(:,datastartmv{i}+k-1)=us(:,k);

    end


    
  

   
    % for l=datastartmv{i}:fullstep
    %     xautofull(:,l)=xtfull(:,l);
    % end
end
xout   = xtfull;%(outdof,:);
% invariance PDE error
% res   = invariance_PDE_residual(DS,tspan,zt,ut);
% ratio = ratio_ext2int(DS,tspan,zt,ut);
%% output
traj      = struct();
traj.time = tspan;
traj.us   = real(utfull);
%traj.etat = etat;
traj.xit  = xisfull;
traj.xt   = real(xout);
traj.xauto=xauto;
%traj.mlus=usmatlab;
%traj.mlxt=real(P*ximatlab);

end


