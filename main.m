%Start script
%%
%modeling Euler-Bernoulli beam in MFE
clear all;
%Parametros
N = 5;%Nelementos
E = 7e10;%modulo de young
Le = 5; %comprimento da barra
le = Le/N; %comprimento do elemento
rho = 2700;%densidade da barra
h = 0.5;%comprimento
b = 0.5;%largura
Ae = b*h;%area
I = (b*h*(b^2 + h^2))/24;%inercia

%Obtendo as matrizez de massa e rigidez
[K3,M3] = beam_ef_creating(N,rho,E,Ae,I,Le);
%amortecimento proporcional
alfa=1e-2;
beta=1e-5;
C3=alfa*M3 + beta*K3;
[m,n]=size(M3);
%modelelando em elementos finitos
%% random field parameters
%% initial
d_KL  = 4;           % number of terms in the expansion
ns = 100;
mu    = 1;           % mean of the field
sigma = 1;         % standard deviation of the field
l_c   = 1/2;           % correlation length x  c = 1/b = 1/2 = 2
% domain of definition
a    = 0.25;
b    = 0.5;
D    = b-a;
n    = 3e2;
xnod = linspace(-1/2,1/2,n);
uo = 0.001;

cov_kernel = @(x1,x2) sigma^2* exp(-abs(x1-x2)/l_c); %fun��o de covari�ncia 
[eigval, eigvec, eigfun, wn,alphaVec] = KL_analytical(d_KL, l_c, sigma, xnod, 'true');%
rhos=rho+5*rho/100.*randn(100,1);% varia��o da densidade, e � randomico
%%

for p = 1: ns %simula��o de monte carlo 
    
    [K_s,M_s] = beam_ef_creating(N,rhos(p),E,Ae,I,Le); % cria a barra em elementos finitos porem levando em conta o par�metro incerto
    M_ie = zeros(m);
    for i = 1:d_KL
        M_ie = M_ie + sqrt(alphaVec(i).*eigval(i))*M_s; %cria a matriz rand�mica
    end
    Mt = M3 + M_ie; %ele soma a matriz rand�mica com a matriz global
    
    
    [U,t]=time_response(Mt,K3,C3,1,uo); %resposta temporal em rela��o a deflex�o da viga
    
    figure(3)
    plot(t,U(:,end-1));
    xlabel('Tempo[s]','Interpreter','latex','fontsize',16)
    ylabel('$u_{n}$[m]','Interpreter','latex','fontsize',16)
    hold on
end





