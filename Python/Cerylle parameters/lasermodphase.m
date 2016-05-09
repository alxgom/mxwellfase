clear all
close all 
clear variables

%global sigma delta lambda beta good%m omega gamma GA good k C GR GRP ktps
global k mu dphim delta gamma do omega dphio vecd indo

indfig=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parametres initiaux 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% parametre des equations 
gammapar=1e8;
k=9e6/gammapar;
omega=0.001*2*pi;  %0.00102 periode en temps renormalisé
mu=0.25e-4;
dphim=0.02;%0.02;
delta=0;%0.000001;%0.01;
gamma=2.5e4/gammapar;
do=2/mu*k;
dphio=0;

% gammapar=1e8;
% k=9e6/gammapar;
% omega=0.0005*2*pi;  %0.00102 periode en temps renormalisé
% mu=0.25e-4;
% dphim=0.00100;
% delta=0;%0.01;
% gamma=0.00027*k%2.5e4/gammapar;
% do=1.3/mu*k;
% dphio=0;
indo=0;

%%%%%%%%%%%%%% parametre pour les temps d'analyse

tempsmin=0;%1e3; % temps a partir duquelle commence l'enregistrement
tempsmax=500e-6*gammapar;%2e5; %5e6*2*pi; % temps final pour le calcul de l'integration

%%%%%%%%%%%%%% parametre pour la variation de la modulation

%%%%%%%%%%%%% parametres pour la tolerance dans la methode de resolution
%%%%%%%%%%%%%% ODE

Rtol=1e-6; % ressere la reponse dans l espace des phases
Atol=1e-46; % tres important pour la reussite de lode


%%%%%%%%%%% Calcul renormalisées

disp('équations renormalisées + temps')

%Iini=25000; %condition initiale sur I
%Nini=122222; % condition initiale pour la phase

indlll=0;

tmin=0;
tmax=tempsmax;
condini=[0.75  0.1 0.75 0.1 2600 0.5 2600 00.5 3600];

%for lambda=35:1:50
%lambda

good=0; %%%%%%%% parametre pour ode45 pour savoir s'il faut modifier Atol
indtest=0;
indlll=indlll+1;

while good==0
    
    options = odeset('RelTol',Rtol,'AbsTol',Atol); %fixe Rtol et Atol en fonction des parametres initiaux
    
    [T1,Y1] = ode45('equa',[0 tmax],condini,options); %lance la resolution
    
    if good==0  %diminue la resolution si cela n a pas convergé
        if Atol >1e-46
            Atol=Atol/100;
            indtest=indtest+1;
            disp('diminution de la resolution');
        else
            good=1;
        end
    end
    if indtest>5 % si 5 diminution de Atol pas suffisant diminution de Rtol
        indtest=0;
        Rtol=Rtol/10;
    end
    
    
end

Inten=(Y1(:,1)+i*Y1(:,2)).*(Y1(:,1)-i*Y1(:,2))+(Y1(:,3)+i*Y1(:,4)).*(Y1(:,3)-i*Y1(:,4));
EX=(Y1(:,1)+i*Y1(:,2)).*(Y1(:,1)-i*Y1(:,2));
EY=(Y1(:,3)+i*Y1(:,4)).*(Y1(:,3)-i*Y1(:,4));

figure(indfig)
hold on
plotyy(T1/gammapar*1e6,Inten,T1/gammapar*1e6,dphio+dphim*cos(omega*T1))
%plot(T1,(Y1(:,1)+i*Y1(:,2)).*(Y1(:,1)-i*Y1(:,2))+(Y1(:,3)+i*Y1(:,4)).*(Y1(:,3)-i*Y1(:,4)),'-b')
xlabel('real time')

figure(indfig+1)
hold on
plot(T1/gammapar*1e6,Inten,'-b')
xlabel('real time')

figure(indfig+2)
hold on
plot(Y1(:,9),Inten,'-b')
xlabel('population')
ylabel('intensity')

figure(indfig+3)
hold on
plot(T1/gammapar*1e6,Y1(:,1),'-+b',T1/gammapar*1e6,Y1(:,3),'-k')
legend('real Ex','real Ey')
xlabel('real time')

figure(indfig+4)
hold on
plot(T1/gammapar*1e6,Y1(:,2),'-+b',T1/gammapar*1e6,Y1(:,4),'-k')
legend('Im Ex','Im Ey')
xlabel('real time')

figure(indfig+5)
hold on
plot(T1,Inten,'-b')
xlabel('renormalized time')

figure(indfig+6)
hold on
plot(T1/gammapar*1e6,Inten,'-b',T1/gammapar*1e6,EX,'-k',T1/gammapar*1e6,EY,'-r')
legend('I','|Ex|','|Ey|')
xlabel('real time')

