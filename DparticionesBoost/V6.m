% ========================================================================
% Fronteras de Estabilidad para convertidor CD-CD Boost
% Autor : Erick Moreno Negrete
% Version : 1.6
% ========================================================================

clc
format longG



% ========================================================================
%                   Parametros del Sistema
% ========================================================================

L = 2.7648e-3;
C = 1.666667e-6;
R = 144;
E = 48;
d = 0.6;
dp = 1-d;
Vo = E/(1-d); % Equilibrio de Voltaje
IL = E/(R*((1-d)^2)); % Equilibrio de Corriente

% ========================================================================
%                     Modelo de pequeña señal
% ========================================================================

A = [   0     -dp/L;
     dp/C  -1/(C*R)];

F = [Vo/L;
     -IL/C];
 
C1 = [1 0];
C2 = [0 1];

D = 0;

% ========================================================================
%                    Funciones de transferencia 
% ========================================================================

numC = [0 Vo/L (2*dp*IL)/(L*C)];
numV = [0 -IL/C (Vo*dp)/(L*C)];
den =  [1 1/(R*C) (dp^2)/(L*C)];

G0 = tf(numV,den);  % Voltaje/d
G1 = tf(numC,den);  % Corriente/d
G2 = tf(numV,numC); % Voltaje/Corriente

% elementos de diseño
a1 = -numV(2);
a2 = numV(3);
b1 = numC(2);
b2 = numC(3);


% ========================================================================
%       Funcion de transferencia lazo interno de corriente 
% ========================================================================

% kp = 0.136597035040431;
% ki = 557.831189710611;

kp = 0.213453038674033;
ki = 1681.51315789474;

nLC = [(kp*Vo)/L ((2*dp*IL*kp)+(Vo*C*ki))/(L*C) (2*dp*IL*ki)/(L*C)];
dLC = [1 ((kp*Vo*R*C)+L)/(R*L*C) ...
       ((2*dp*IL*kp)+(dp^2)+(Vo*C*ki))/(L*C) ...
       (2*dp*IL*ki)/(L*C)];
   
GCcl = tf(nLC,dLC);

% elementos de diseño
n1 = nLC(1);
n2 = nLC(2);
n3 = nLC(3);

d1 = dLC(2);
d2 = dLC(3);
d3 = dLC(4);


% ========================================================================
%       Funcion de transferencia lazo externo de voltaje (lazo abierto)
% ========================================================================

nVLA = [-n1*a1 (n1*a2)-(n2*a1) (n2*a2)-(n3*a1) n3*a2];
dVLA = [b1 b2+(d1*b1) (d1*b2)+(d2*b1) (d2*b2)+(d3*b1) d3*b2];

% elementos de diseño
nV1 = nVLA(1);
nV2 = nVLA(2);
nV3 = nVLA(3);
nV4 = nVLA(4);

dV1 = dVLA(1);
dV2 = dVLA(2);
dV3 = dVLA(3);
dV4 = dVLA(4);
dV5 = dVLA(5);

% ========================================================================
%     Caracterizacion de las fronteras de estabilidad lazo externo
% ========================================================================

% Fronteras de sigma

w = linspace(0,17000,10000);
sb = 1i.*w;

s1 = -500+(1i.*w);
s1n = -500;
    
s2 = -1000+(1i.*w);
s2n = -1000;
    
s3 = -4000+(1i.*w);
s3n = -4000;

%Polinomios en s=jw; fronteras de estabilidad
Nds0 = nV1*sb.^3 + nV2*sb.^2 + nV3*sb + nV4;
Dds0 = dV1*sb.^4 + dV2*sb.^3 + dV3*sb.^2 + dV4*sb + dV5;

%Polinomios en s=-100+jw; 
Nds1 = nV1*s1.^3 + nV2*s1.^2 + nV3*s1 + nV4;
Dds1 = dV1*s1.^4 + dV2*s1.^3 + dV3*s1.^2 + dV4*s1 + dV5;

%Polinomios en s=-100; 
Nds1n = nV1*s1n.^3 + nV2*s1n.^2 + nV3*s1n + nV4;
Dds1n = dV1*s1n.^4 + dV2*s1n.^3 + dV3*s1n.^2 + dV4*s1n + dV5;

%Polinomios en s=-500+jw; 
Nds2 = nV1*s2.^3 + nV2*s2.^2 + nV3*s2 + nV4;
Dds2 = dV1*s2.^4 + dV2*s2.^3 + dV3*s2.^2 + dV4*s2 + dV5;
    
%Polinomios en s=-500; 
Nds2n = nV1*s2n.^3 + nV2*s2n.^2 + nV3*s2n + nV4;
Dds2n = dV1*s2n.^4 + dV2*s2n.^3 + dV3*s2n.^2 + dV4*s2n + dV5;

%Polinomios en s=-1000+jw; 
Nds3 = nV1*s3.^3 + nV2*s3.^2 + nV3*s3 + nV4;
Dds3 = dV1*s3.^4 + dV2*s3.^3 + dV3*s3.^2 + dV4*s3 + dV5;
    
%Polinomios en s=-1000; 
Nds3n = nV1*s3n.^3 + nV2*s3n.^2 + nV3*s3n + nV4;
Dds3n = dV1*s3n.^4 + dV2*s3n.^3 + dV3*s3n.^2 + dV4*s3n + dV5;


vKp0 = -real(Dds0./Nds0);
vKi0 = w.*imag(Dds0./Nds0);
rKp0 = linspace(-0.05,0.03,1000);
rKi0 = linspace(0,0,1000);
rKi0y = linspace(-1,320,1000);

figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'GridLineStyle','--')

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#2f3e46';
color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

% sigma = 0
plot(vKp0,vKi0,'Color',color1)
hold on
plot(rKp0,rKi0,'b--','HandleVisibility','off')
plot(rKi0,rKi0y,'b--','HandleVisibility','off')

% sigma = 100
vKp1 = -real(Dds1./Nds1);
vKi1 = w.*imag(Dds1./Nds1);

%recta 
rKp1 = linspace(-0.05,0.03,1000);
rKi1 = (50*rKp1)+50*(Dds1n/Nds1n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#354f52';
color2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp1,vKi1,'--','Color',color2)
plot(rKp1,rKi1,'--','Color',color2,'HandleVisibility','off')

% sigma = 500
vKp2 = -real(Dds2./Nds2);
vKi2 = w.*imag(Dds2./Nds2);
        
%recta 
rKp2 = linspace(-0.05,0.03,1000);
rKi2 = (1000*rKp2)+1000*(Dds2n/Nds2n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#52796f';
color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp2,vKi2,':','Color',color3)
plot(rKp2,rKi2,':','Color',color3,'HandleVisibility','off')

% sigma = 1000
vKp3 = -real(Dds3./Nds3);
vKi3 = w.*imag(Dds3./Nds3);
        
%recta 
rKp3 = linspace(-0.05,0.03,1000);
rKi3 = (1000*rKp3)+1000*(Dds3n/Nds3n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#84a98c';
color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp3,vKi3,'-.','Color',color3)
plot(rKp3,rKi3,'-.','Color',color3,'HandleVisibility','off')

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#f72585';
color4 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(0.5e-3,50,'p','Color',color4)

legend('$\sigma=0$','$\sigma<100$','$\sigma<500$','$\sigma<1000$',...
    '$(k_p^v,k_i^v)$', 'interpreter', 'latex','Location', 'Best')
axis([-0.04 0.04 0 250]);
xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');


