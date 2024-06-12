function [S,Svoigt,Pvoigt,sigma,epsilon,strainEnergy,matStiffDensity] = elastic_neoHookean(F,s,materialCfg)

DOF=length(F);
if DOF==2
    Voigt=[1,1;2,2;1,2];
    len=3;
else
    Voigt=[1,1;2,2;3,3;2,3;1,3;1,2];
    len=6;
end
C=F'*F;
[v,a]=eig(C);
a=diag(a);
[a,I]=sort(a);
v=v(:,I);
J=prod(sqrt(a));
E=1/2*(C-eye(DOF,DOF));


W1=0;
mu=zeros(1,DOF);
for i=1:DOF
    Pa=a(i);
    mu(i)=MUsign(a(i),materialCfg);
    W1=W1+mu(i)/2*(Pa-1-log(Pa));
end

K=Ksign(J,materialCfg);

W1=W1+K/2*log(J)^2;

strainEnergy(1)=W1; % positive part
strainEnergy(2)=0; % negative part
S=PK2(J,K,mu,a,v,C);
S=s*S;
Kp=Kmat(J,K,mu,a,v,C,0);
matStiffDensity=s*Kp;


Svoigt=zeros(len,1);
for a=1:len
    i=Voigt(a,1);
    j=Voigt(a,2);
    Svoigt(a)=S(i,j);
end

P=S*F';
Pvoigt=zeros(len,1);
for a=1:len
    i=Voigt(a,1);
    j=Voigt(a,2);
    Pvoigt(a)=P(i,j);
end

epsilon=zeros(len,1);
for a=1:len
    i=Voigt(a,1);
    j=Voigt(a,2);
    epsilon(a)=E(i,j);
    epsilon(DOF+1:end)=2*epsilon(DOF+1:end);
end

T=F*S*F'/J;
sigma=zeros(len,1);
for a=1:len
    i=Voigt(a,1);
    j=Voigt(a,2);
    sigma(a)=T(i,j);
end
end

function y=Psign(x)
if x>1
    y=x;
else
    y=1;
end
end

function y=Nsign(x)
if x<=1
    y=x;
else
    y=1;
end
end

function k=Ksign(J,materialCfg)
if J>1
    k=materialCfg.K_pos;
else
    k=materialCfg.K_neg;
end
end

function mu=MUsign(lambda,materialCfg)
if lambda>1
    mu=materialCfg.mu_pos;
else
    mu=materialCfg.mu_neg;
end
end

function S=PK2(J,K,mu,a,v,C)
DOF=length(C);
S=zeros(DOF,DOF);
for k=1:DOF
    S=S+mu(k)*(1-1/a(k))*kron(v(:,k),v(:,k)');
end
S=S+K*log(J)*C^(-1);
end

function D=Kmat(J,K,mu,a,v,C,sig)
DOF=length(C);
invC=inv(C);
if DOF==2
    Voigt=[1,1;2,2;1,2];
    len=3;
else
    Voigt=[1,1;2,2;3,3;2,3;1,3;1,2];
    len=6;
end
% https://mathoverflow.net/questions/229425/derivative-of-eigenvectors-of-a-matrix-with-respect-to-its-components
% Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill’s family of generalized strain tensors
D=zeros(len,len);
for p=1:len
    for q=1:len
        if p>=q
            i=Voigt(p,1);
            j=Voigt(p,2);
            m=Voigt(q,1);
            n=Voigt(q,2);
            for k=1:DOF
                pVikpCmn=0;
                pVjkpCmn=0;
                for b=1:DOF
                    if b~=k && ~isinf(1/(a(k)-a(b))) && a(k)~=1
                        pVikpCmn=pVikpCmn+1/2/(a(k)-a(b))*v(i,b)*(v(m,k)*v(n,b)+v(m,b)*v(n,k));
                        pVjkpCmn=pVjkpCmn+1/2/(a(k)-a(b))*v(j,b)*(v(m,k)*v(n,b)+v(m,b)*v(n,k));
                    end
                end
                D(p,q)=D(p,q)+2*mu(k)*(1/a(k)^2*v(m,k)*v(n,k)*v(i,k)*v(j,k)+(1-1/a(k))*(pVikpCmn*v(j,k)+v(i,k)*pVjkpCmn))*myHeaviside((a(k)-1),sig);
            end
            D(p,q)=D(p,q)+(K*invC(m,n)*invC(i,j)-K*log(J)*(invC(i,m)*invC(j,n)+invC(i,n)*invC(j,m)))*myHeaviside((J-1),sig);
        end
    end
end
for p=1:len
    for q=1:len
        if p<q
            D(p,q)=D(q,p);
        end
    end
end
end

function y=myHeaviside(x,sig)
if sig==0
    y=1;
elseif sig>0
    if x==0
        y=1;
    else
        y=heaviside(x);
    end
elseif sig<0
    if x==0
        y=0;
    else
        y=heaviside(-x);
    end
end
end