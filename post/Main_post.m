load data1.mat
imax        = Nx;
omega       = NPSE(4);
beta        = NPSE(5);
Re0         = parameter(7);
Ma          = parameter(3);
Te          = parameter(4);



zz=sqrt(-1);

x=X;
R=sqrt(x*Re0);

alf1 = reshape(alf,[2*n_max+1,m_max+1,imax]);
fai  = reshape(Fai,[5*N,2*n_max+1,m_max+1,imax]);

JF=zeros(imax,m_max+1,2*n_max+1);
for k=1:2*n_max+1    
for j=1:m_max+1
for i=1:imax
        if i>=2
            JF(i,j,k)=trapz(x(1:i),alf1(k,j,1:i));
        else
            JF(i,j,k)=0;
        end
end
end
end


rho=zeros(Ny,imax,m_max+1,2*n_max+1);
U=zeros(Ny,imax,m_max+1,2*n_max+1);
V=zeros(Ny,imax,m_max+1,2*n_max+1);
W=zeros(Ny,imax,m_max+1,2*n_max+1);
T=zeros(Ny,imax,m_max+1,2*n_max+1);
for k=1:2*n_max+1
    for j=1:m_max+1
        for i=1:imax
            for l=1:Ny
            rho(l,i,j,k)=fai(5*(l-1)+1,k,j,i)*exp(zz*JF(i,j,k));
            U(l,i,j,k)=fai(5*(l-1)+2,k,j,i)*exp(zz*JF(i,j,k));
            V(l,i,j,k)=fai(5*(l-1)+3,k,j,i)*exp(zz*JF(i,j,k));
            W(l,i,j,k)=fai(5*(l-1)+4,k,j,i)*exp(zz*JF(i,j,k));
            T(l,i,j,k)=fai(5*(l-1)+5,k,j,i)*exp(zz*JF(i,j,k));
            end
        end
    end
end

%figure 10
umax1=max(abs(U(:,:,2,1)))*sqrt(2);
umax2=max(abs(U(:,:,3,1)))*sqrt(2);

plot(R,umax1,'k',R,umax2,'k.')
hold on

plot(R,umax1,'b',R,umax2,'b.')
hold on

%d1=csvread('0.25-F.csv') ; d2=csvread('0.25-2F.csv') ;
%scatter(d1(:,1),d1(:,2),20,'filled','r');hold on
%scatter(d2(:,1),d2(:,2),20,'filled','b');hold on
%legend('0.25%-F','0.25%-2F','0.25%-F-文献','0.25%-2F-文献')
%xlabel('R')
%ylabel('u^{,}_{max}')

%d1=csvread('0.3-F.csv') ; d2=csvread('0.3-2F.csv') ;
%scatter(d1(:,1),d1(:,2),20,'filled','b');hold on
%scatter(d2(:,1),d2(:,2),20,'filled','b');hold on
%legend('0.30%-F','0.30%-2F','0.30%-F-文献','0.30%-2F-文献')
%xlabel('R')
%ylabel('u^{,}_{max}')


legend('0.25%-F','0.25%-2F','0.30%-F','0.30%-2F')
xlabel('R')
ylabel('u^{,}_{max}')


%figure 11

eta=y*400/793.9773296511684;
d11=csvread('0F.csv') ; d22=csvread('1F.csv') ; d33=csvread('2F.csv') ; d44=csvread('3F.csv'); 

%1
plot(U(:,57,1,1),eta,'k')
hold on
axis([-0.0015 0.0015 0 8])
scatter(d11(:,1),d11(:,2),20,'filled','b');
legend('U-0F','U-0F-文献')
ylabel('\eta')
xlabel('u')

%2
plot(abs(U(:,57,2,1))*sqrt(2),eta,'k')
hold on
axis([0 0.02 0 8])
scatter(d22(:,1),d22(:,2),20,'filled','b');
legend('U-1F','U-1F-文献')
ylabel('\eta')
xlabel('u^{,}')


%3
plot(abs(U(:,57,3,1))*sqrt(2),eta,'k')
hold on
axis([0 0.002 0 8])
scatter(d33(:,1),d33(:,2),20,'filled','b');
legend('U-2F','U-2F-文献')
ylabel('\eta')
xlabel('u^{,}')

%4
plot(abs(U(:,57,4,1))*sqrt(2),eta,'k')
hold on
axis([0 0.00025 0 8])
scatter(d44(:,1),d44(:,2),20,'filled','b');
legend('U-3F','U-3F-文献')
ylabel('\eta')
xlabel('u^{,}')
