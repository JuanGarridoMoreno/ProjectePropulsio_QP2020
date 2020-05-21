% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

% El codi determina l'impuls i empenta espec�fics a partir del set de
% par�metres inicial

clear;
close all;
clc;

% Introducci� de les dades del motor a analitzar. Totes elles expressades
% segons el sistema internacional d'unitats.
M_0=0.85; z=11000; T_SL=288.15; P_SL=101325;
pi_d=0.98; eta_f=0.89; eta_lpc=0.88; eta_hpc=0.86;
pi_b=0.99; eta_b=0.99; eta_hpt=0.91; eta_lpt=0.92;
eta_mH=0.993; eta_mL=0.997; pi_np=0.99; pi_ns=0.99;
T_t4=1450; gamma_c=1.4; gamma_t=1.3; R_g=287; h=43000000; g=9.81;

%% C�LCULS PREVIS

% Es calculen els coeficients de pressi� per l'aire del compressor (Cp_c) i
% de la turbina (Cp_t), segons les gammes depenent de si ha passat o no per
% la cambra de combusti�.
Cp_t=(gamma_t/(gamma_t-1))*R_g;
Cp_c=(gamma_c/(gamma_c-1))*R_g;

% S'obtenen tamb� els valors de pressi� i temperatura del flux lliure a
% l'altura de funcionament. Amb aquests es calula la velocitat d'entrada de
% l'aire al motor.
[T_0,P_0]=atmosferaISA(z,T_SL,P_SL);
u_0=M_0*sqrt(gamma_c*R_g*T_0);

%% Set de par�metres

% Es proposa un grup de par�metres que d�na relativament bons resultats per
% l'impuls i l'empenta espec�fics.
alpha=8;
pi_f=1.6;
pi_lpc=4.3;
pi_hpc=7.5;
OPR=pi_lpc*pi_hpc;

% Per tal de fer l'an�lisi del motor, es divideix el proc�s en les etapes
%que es troba l'aire al seu pas aquest. La notaci� emprada correspon a
%la que apareix a l'esquema del motor.

%% (0)-(2): INLET

theta_0=1+0.5*(gamma_c-1)*M_0^2; % theta_0=T_t0/T_0
delta_0=theta_0^(gamma_c/(gamma_c-1)); % delta_0=P_t0/P_0

% Es calcula la temperatura d'estancament del flux lliure, que roman
% constant durant la presa d'aire.

T_t0=T_0*theta_0;
T_t2=T_t0;

% Es calculen tamb� les pressions d'estancament. Aquesta �s menor a la presa
% d'aire degut a l'efici�ncia d'aquesta.

P_t0=P_0*delta_0;
P_t2=P_t0*pi_d;

%% (2)-(2.5): LPC

% Es calcula la pressi� al compressor de baixa pressi� gr�cies al rendiment
% isentr�pic d'aquest. A partir d'aquest, i de la gamma de l'aire s'obt�
% l'efici�ncia isentr�pica que permet con�ixer la temperatura.

P_t25=P_t2*pi_lpc;
tau_lpc=1+(pi_lpc^((gamma_c-1)/gamma_c)-1)/eta_lpc;
T_t25=T_t2*tau_lpc;

%% (2.5)-(3): HPC

% El procediment de c�lcul al compressor d'alta pressi� �s el mateix que al
% de baixa, partint de la temperatura i la pressi� d'aquest �ltim.

P_t3=P_t25*pi_hpc;
tau_hpc=1+(pi_hpc^((gamma_c-1)/gamma_c)-1)/eta_hpc;
T_t3=T_t25*tau_hpc;

%% (3)-(4): Cambra de Combusti�

% Amb les propietats f�siques de l'aire a l'entrada i la sortida d'aquesta i
% l'efici�ncia isentr�pica es pot con�ixer la relaci� de cabal de
% combustible respecte el cabal d'aire. I amb la relaci� de pressions
% d'estancament, la pressi� a la sortida.

f=(Cp_t*T_t4-Cp_c*T_t3)/(eta_b*h-Cp_t*T_t4);
P_t4=P_t3*pi_b;

%% (4)-(4.5): HPT

% Per obtenir la relaci� de temperatures a la turbina d'alta pressi� s'usen
% l'efici�ncia mec�nica en aquesta, aix� com les temperatures al compressor
%i a l'entrada de la HPT, tenint en compte la pres�ncia de combustible. A
% partir d'aquesta relaci�, ja es poden calcular la relaci� de pressions i
% la temperatura i pressi� a la sortida d'aquest element.

tau_hpt=1-(1/eta_mH)*(Cp_c/Cp_t)*(1/(1+f))*(T_t25/T_t4)*(tau_hpc-1);
pi_hpt=(1-(1-tau_hpt)/eta_hpt)^(gamma_t/(gamma_t-1));
T_t45=T_t4*tau_hpt;
P_t45=P_t4*pi_hpt;

%% (4.5)-(5): LPT

%La turbina de baixa pressi� segueix el mateix concepte i procediment que
%la turbina d'alta pressi�.

tau_f=1+(pi_f^((gamma_c-1)/gamma_c)-1)/eta_f;
tau_lpt=1-(1/eta_mL)*(Cp_c/Cp_t)*(1/(1+f))*(T_t0/T_t4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_f-1));
pi_lpt=(1-(1-tau_lpt)/eta_lpt)^(gamma_t/(gamma_t-1));
T_t5=T_t45*tau_lpt;
P_t5=P_t45*pi_lpt;

%% (5)-(9): PRIMARY NOZZLE

%Les temperatures d'estancament de la tovera prim�ria es mantenen
%constants, mentre que la pressi� disminueix lleugerament a causa de
%l'efici�ncia de la tovera.

P_t9=P_t5*pi_np;
T_t9=T_t5;

%Relaci� de pressi� d'estancament a la tovera i a l'entrada del motor.
pi_primari_total=P_t9/P_0; 

% Es realitza la hip�tesi de tovera adaptada
M_9=sqrt(2/(gamma_t-1)*((pi_primari_total)^((gamma_t-1)/gamma_t)-1));

if M_9<1 % Hip�tesi certa
    P_9=P_0;
    pi_primari=1;
else     %Hip�tesi falsa: tovera ofegada, el Mach no pot ser superior a la unitat.
    M_9=1;
    pi_primari=pi_primari_total*(2/(gamma_t+1))^(gamma_t/(gamma_t-1)); 
    % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que �s fluxe secundari
    P_9=pi_primari*P_0;
end
T_9=T_t9/(1+0.5*(gamma_t-1)*M_9^2);
u_9=M_9*sqrt(gamma_t*R_g*T_9);

%% (0)-(1.3): SECONDARY INLET

% A la presa second�ria d'aire, la temperatura d'estancament �s major que la
% del flux lliure.
T_t13=T_0*theta_0*tau_f;

%% (1.3)-(1.9): SECONDARY NOZZLE

% Com l'aire no experimenta cap altre proc�s despr�s de la presa second�ria,
% la temperatura d'estancament a la tovera second�ria �s la mateixa que a la
% presa. Tot i aix� cal diferenciar si aquesta es troba xocada o no.
T_t19=T_t13;

% pi_secundari_total=P_t19/P_0
pi_secundari_total=delta_0*pi_d*pi_f*pi_ns; 

% Es realitza la hip�tesi de tovera adaptada i es calcula el Mach
M_19=sqrt(2/(gamma_c-1)*((pi_secundari_total)^((gamma_c-1)/gamma_c)-1));

if M_19<1 % Hip�tesi certa
    P_19=P_0;
    pi_secundari=1;
else      % Hip�tesi falsa: tovera ofegada, el Mach no pot ser superior a la unitat.
    M_19=1;
    % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que �s fluxe secundari
    pi_secundari=pi_secundari_total*(2/(gamma_c+1))^(gamma_c/(gamma_c-1)); 
    P_19=pi_secundari*P_0;
end
T_19=T_t19/(1+(gamma_c-1)/2*M_19^2);
u_19=M_19*sqrt(gamma_c*R_g*T_19);

%% C�lcul de l'empenyiment espec�fic

% Es plantegen les equacions d'empenyiment dels dos tipus de flux i es fan
% relatives al flux d'aire que veu el motor. L'empenta total que veu el
% motor �s la suma de les dues.

Specific_Thrust_primari=(1+f)*u_9-u_0+(1+f)*T_9/u_9*R_g*(1-1/pi_primari);
Specific_Thrust_secundari=alpha*(u_19-u_0)+alpha*(T_19*R_g)/u_19*(1-1/pi_secundari);
Specific_Thrust=Specific_Thrust_primari+Specific_Thrust_secundari;

%% C�lcul de l'impuls espec�fic

% L'impuls espec�fic es calcula com la relaci� entre la for�a produ�da i el
% combustible emprat.
I_sp=Specific_Thrust/(f*g);
if imag(I_sp)~=0||imag(Specific_Thrust)~=0
    Specific_Thrust=0;
    I_sp=0;
end

% Impressi� de resultats
I_sp
Specific_Thrust

P_0;
P_9;
P_19;