function [Admv,Bdmv,Uupmv,datastartmv,datalengthmv,DMDtraj]=DMD1(datalength,datastart,r,windownum,X,urand,reduceornot)
close all
%

%tcdim=100;
% datalength=9000;
% datastart=1;
% %p=tcdim+2;
% r=20;
% windownum=2


Admv=cell(1,windownum);
Bdmv=cell(1,windownum);
Uupmv=cell(1,windownum);
datastartmv=cell(1,windownum);
datalengthmv=cell(1,windownum);
%%%%%%%%%%%ron
if reduceornot==0
    r=size(X,1);
end
xtrkre=[];
yt=[];
%X=xfullwRc(1:end-1,:)'-yRec{2};

%p=45;

thresmax=3;


for i=1:windownum


    Xhalf=X(:,datastart:datastart+datalength);
    %X=trajs3{1}.zt(1:tcdim,:);
    u=urand(:,datastart:datastart+datalength-1);
    N=size(Xhalf,2)-1;
    maxdata=thresmax*max(abs(Xhalf(1,:)));



    X1 = Xhalf(:,1:end-1);
    X2 = Xhalf(:,2:end);
    Gama=u(:,1:size(X1,2));
    Omega=[X1;Gama];

    [U,S,V] = svd(Omega,'econ');
    p=size(Omega,1);

    Utilde = U(:,1:p);
    Stilde = S(1:p,1:p);
    Vtilde = V(:,1:p);

    U1tilde=Utilde(1:size(X1,1),:);
    U2tilde=Utilde(size(X1,1)+1:size(Utilde,1), :);

    A_=X2*Vtilde*pinv(Stilde)*(U1tilde.');
    B_=X2*Vtilde*pinv(Stilde)*(U2tilde.');%3.17

    [Uup,Sup,Vup]=svd(X2,'econ');

    Uuptd = Uup(:,1:r);
    Suptd = Sup(1:r,1:r);
    Vuptd = Vup(:,1:r);

    Atilde=(Uuptd.')*A_*Uuptd;%3.18
    Btilde=(Uuptd.')*B_;

    % [W,eigs] = eig(Atilde);
    % phi=X2*Vtilde*inv(Stilde)*(U1tilde.')*Uuptd*W; %3.21
   

    Ad{i}=A_;
    Bd{i}=B_;
        
    Admv{i}=Atilde;
    Bdmv{i}=Btilde;
    Uupmv{i}=Uuptd;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%reduced or not
    if reduceornot==0
    Admv{i}=A_;
    Bdmv{i}=B_;
    Uupmv{i}=eye(size(X,1),size(X,1));
    end

    xt=zeros(size(Ad{i},2),N);
    %yt=[];

    %ut=zeros(2,N);
    ut=zeros(size(urand,1),N);


    epsino=zeros(r,N);
    ytrk=zeros(1,N);
    utrk=zeros(size(urand,1),N);

    X0=X(:,datastart);
    xt(:,1)=X0;

    epsino(:,1)=pinv(Uupmv{i})*X0;
   
     for k =1:datalength-1
        ut(:,k)=u(:,k);
        utrk(:,k)=u(:,k);

        xt(:,k+1)=Ad{i}*xt(:,k)+Bd{i}*ut(:,k);
        epsino(:,k+1)=Admv{i}*epsino(:,k)+Bdmv{i}*ut(:,k);
        %epsino(:,k+1)=Admv{i}*epsino(:,k)+1*ut(:,k);

     end

    xtrkrenew=Uupmv{i}*epsino;
    
    xtrkre=[xtrkre,xtrkrenew];
    

    datalengthmv{i}=datalength;
    datastartmv{i}=datastart;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % datalength1=0; 
    % datastart1=datastart+datalength;
    % retnum=0
    % while max(abs(xtrkrenew(1,:)))>maxdata
    % 
    %     retnum=retnum+1;
    % 
    %     itepointlength=100;
    %     if retnum==1
    %         datalength1=datalength
    %     end
    % 
    %     datalength1=datalength1+itepointlength;
    %     if datastart1>datalength1  %不能减到负数
    %     datastart1=datastart1-itepointlength;
    %     else
    %         datalength1=datalength1-itepointlength;
    %         break
    %     end
    % 
    %     Xhalf1=X(:,datastart1:datastart1+datalength1);
    %     u=0.1*traj_closed.ut(:,datastart1:datastart1+datalength1-1);
    %     N=size(Xhalf1,2)-1;
    % 
    % 
    % 
    % 
    %     X1 = Xhalf1(:,1:end-1);
    %     X2 = Xhalf1(:,2:end);
    %     Gama=u(:,1:end);
    %     Omega=[X1;Gama];
    % 
    %     [U,S,V] = svd(Omega,'econ');
    % 
    %     Utilde = U(:,1:p);
    %     Stilde = S(1:p,1:p);
    %     Vtilde = V(:,1:p);
    % 
    %     U1tilde=Utilde(1:size(X1,1),:);
    %     U2tilde=Utilde(size(X1,1)+1:size(Utilde,1), :);
    % 
    %     A_=X2*Vtilde*pinv(Stilde)*(U1tilde.')
    %     B_=X2*Vtilde*pinv(Stilde)*(U2tilde.')%3.17
    % 
    %     [Uup,Sup,Vup]=svd(X2,'econ');
    % 
    %     Uuptd = Uup(:,1:r);
    %     Suptd = Sup(1:r,1:r);
    %     Vuptd = Vup(:,1:r);
    % 
    %     Atilde=(Uuptd.')*A_*Uuptd%3.18
    %     Btilde=(Uuptd.')*B_
    % 
    %     Uupmv{i}=ones(20,20);%Uuptd;
    % 
    % 
    %     Ad{i}=A_;
    %     Bd{i}=B_;
    %     % Admv{i}=Atilde;
    %     % Bdmv{i}=Btilde;
    %     Admv{i}=A_;
    %     Bdmv{i}=B_;
    % 
    % 
    %     xt=zeros(size(Ad{i},2),N);
    %     %yt=[];
    %     ut=zeros(2,N);
    % 
    %     epsino=zeros(r,N);
    %     ytrk=zeros(1,N);
    %     utrk=zeros(2,N);
    % 
    %     X0=X(:,datastart1);
    %     xt(:,1)=X0;
    %     %epsino(:,1)=pinv(Uuptd)*X0;
    %     epsino(:,1)=X0;
    % 
    %     for k =1:datalength1-1
    % 
    %         utrk(:,k)=u(:,k);
    % 
    % 
    %         epsino(:,k+1)=Admv{i}*epsino(:,k)+Bdmv{i}*utrk(:,k);
    % 
    %     end
    %       %xtrkrenew=Uuptd*epsino;
    %       xtrkrenew=epsino;
    %       maxdata=thresmax*max(abs(Xhalf1(1,:)));
    %     if retnum>i || retnum>40
    %         disp("mamba out")
    %         break;
    % 
    %     end
    % 
    % end
    % 
    % 
    % if retnum~=0
    % xtrkre(:,end-datalength1+1:end)=[];
    % xtrkre=[xtrkre,xtrkrenew];
    % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % datalength1=0; 
    % datastart1=datastart;
    % retnum=0
    % while max(abs(xtrkrenew(1,:)))>maxdata
    % 
    %     retnum=retnum+1;
    % 
    %     itepointlength=500;
    %     if retnum==1
    %         datalength1=datalength
    %     end
    % 
    %     datalength1=datalength1+itepointlength;
    %     if datastart1>datalength1  %不能减到负数
    %     datastart1=datastart1-itepointlength;
    %     else
    %         datalength1=datalength1-itepointlength;
    %         break
    %     end
    % 
    %     Xhalf1=X(:,datastart1:datastart1+datalength1);
    %     u=0.1*traj_closed.ut(:,datastart1:datastart1+datalength1-1);
    %     N=size(Xhalf1,2)-1;
    % 
    % 
    %     Admv{i}=Admv{i-1};
    %     Bdmv{i}=Bdmv{i-1};
    %     Uupmv{i}=Uupmv{i-1};
    % 
    % 
    %     xt=zeros(size(Ad{i},2),N);
    %     %yt=[];
    %     ut=zeros(2,N);
    % 
    %     epsino=zeros(r,N);
    %     ytrk=zeros(1,N);
    %     utrk=zeros(2,N);
    % 
    %     X0=X(:,datastart1);
    %     xt(:,1)=X0;
    %     epsino(:,1)=pinv(Uupmv{i-1})*X0;
    %     for k =1:datalength1-1
    % 
    %         utrk(:,k)=u(:,k);
    % 
    % 
    %         epsino(:,k+1)=Admv{i-1}*epsino(:,k)+Bdmv{i-1}*utrk(:,k);
    % 
    %     end
    %       xtrkrenew=Uupmv{i}*epsino;
    %       maxdata=10*max(abs(Xhalf1(1,:)));
    %     if retnum>i || retnum>1
    %         disp("mamba out")
    %         break;
    % 
    %     end
    % 
    % end
    % 
    % datalengthmv{i}=datalength1;
    % datastartmv{i}=datastart1;
    % 
    % if i~=1
    % end 
    % 
    % if retnum~=0
    % xtrkre(:,end-datalength1+1:end)=[];
    % xtrkre=[xtrkre,xtrkrenew];
    % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %xtrkre=[xtrkre,xtrkrenew];
    %yt=[yt,xtrkre];
    datastart=datastart+datalength;

end


ytrkre=xtrkre(1,:);
y=X(1,1:datalength*i);
yt=xtrkre(1,:);
timei2=linspace(1,datalength*i,datalength*i);
timei2=timei2/500;

figure();
plot(timei2,y,'k-','DisplayName','train origin','Color',[0,0,1]);
hold on;
plot(timei2,yt,'r--','DisplayName','train DMD','Color',[1,0,0]);

 %ylim([-10 10])
DMDtraj=cell(1,2);
DMDtraj{1}=timei2;
DMDtraj{2}=yt;
% xlim([0 18000])
% ylim([-5 5])
%plot(timei2,ytrkre,'b','DisplayName','origin','Color',[0,1,1]);

% figure(50)
% time2=linspace(1,500,500)
% plot(trajs3{1}.zt(1,:),'r')
% plot(trajs3{1}.zauto(1,:),'b')
% plot(trajs3{1}.ut(1,:),'g')
end