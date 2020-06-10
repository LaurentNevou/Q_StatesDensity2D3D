%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% last update 4June 2020, lnev %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paul Harrisson
% Quantum Wells, Wires and Dots.
% 4th edition (2016),
% chap 2 : "Solutions to Schrodinger's equation"
% 2.42: "Two-dimensional systems" page 31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Density of states of a two-dimensional electron gas including nonparabolicity"
% J. A. López-Villanueva, F. Gámiz, I. Melchor, and J. A. Jiménez-Tejada
% Journal of Applied Physics 75, 4267 (1994); doi: 10.1063/1.355967
% http://dx.doi.org/10.1063/1.355967
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "On the Thermodynamics of a Two-Dimensional Electron Gas with Non-Parabolic Dispersion"
% G. Gulyamov, B. T. Abdulazizov
% World Journal of Condensed Matter Physics, 2016, 6, 294-299
% DOI: 10.4236/wjcmp.2016.64028
% https://www.scirp.org/pdf/WJCMP_2016111716503321.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Effect of Temperature and Band Nonparabolicity on Density of States of Two Dimensional Electron Gas"
% G. Gulyamov, P. J. Baymatov, B. T. Abdulazizov
% Journal of Applied Mathematics and Physics, 2016, 4, 272-278
% DOI: 10.4236/jamp.2016.42034
% https://www.scirp.org/pdf/JAMP_2016022514420230.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Material parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L     = 20e-9;                       %% quantum well tickness [m]
dz    = 1e-10;                       %% z resolution          [m]
meff0 = 0.067;                       %% effective electron mass
alpha = 0.7;                         %% non-parabolicity parameter [eV-1]
n     = 20;                          %% number of quantum states (MUST be < L/dz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h   = 6.62606896E-34;               %% Planck constant [J.s]
hbar= h/(2*pi);
e   = 1.602176487E-19;              %% electron charge [C]
me  = 9.10938188E-31;               %% electron mass   [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 0:dz:L;                         %% distance vector  [m]
V0 = z*0;                           %% potential vector [eV]
En=linspace( 0 , 0.5, 1000 );         %% Energy vector    [eV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meff3D=meff0*(1+alpha*En).^(1/3) .*(1+2*alpha*En).^(2/3);
ro3D = (1/(2*pi^2)) * ( (2*e*meff3D*me/(hbar^2)).^(3/2) ) .* sqrt( En );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0(1)=max(En);
V0(end)=max(En);
meff=meff0*(1+alpha*En);

if alpha==0
    [E,psi] = Schroed1D_FEM_f(z,V0,meff0,n);
else
%     Eg=1/alpha;
%     EP=Eg/meff0;
%     Dso=0;
%     [E,psi] = Schrod_2bands_Kane_f(z,V0,Eg,EP,Dso,n,0,0,0,0,0);
    
    dE=1e-2;precision=1e-5;
    meff_mat=repmat(meff',[1 length(z)]);
    [E,psi] = Schrod_Nbands_shoot_f(z,V0,meff_mat,n,En,dE,precision);
    
    n=length(E);
end

meff2D=meff0*(1+2*alpha*En);

for i=1:n
    
    ro2D(En>E(i),i)  = e*meff2D(En>E(i))*me/(pi*hbar^2);
    ro2D( En<E(i),i) = 0;
    
end

roo2D=sum(ro2D,2)/L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[10 50 1200 800],'color','w');
subplot(1,1,1,'fontsize',15)
hold on; grid on;

plot(En,ro3D*1e-6,'b','linewidth',2)
plot(En,roo2D*1e-6,'r','linewidth',2)

xlabel('Energy (eV)')
ylabel('Density of states (cm-3.eV-1)')
legend('\color{blue}Bulk (3D)','\color{red}Quantum Well (2D)','location','northwest')

tt=strcat('meff=',num2str(meff0),'; \alpha=',num2str(alpha),'eV-1; L-QW=',num2str(L*1e9),'nm');
title(tt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

ScF=0.05;
for i=1:n
    psi(:,i)=abs(psi(:,i)).^2/max(abs(psi(:,i)).^2)*ScF + E(i); % normalisation for the plotting
end
          
figure('position',[1000 50 900 800],'color','w');
subplot(1,1,1,'fontsize',15)
hold on; grid on;

for i=1:n
    plot(z*1e9,psi(:,i),'r')
end

xlabel('z (nm)')
ylabel('Energy (eV)')
