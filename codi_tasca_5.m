% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

%A partir dels paràmetres de l'apartat 3, s'estudia les prestacions del
%turbofan en funció de la relació de compressions del fan. En tot moment,
%s'imposa tovera adaptada, tant en la primària com en la secundària.

clear
close all
clc

% Introducció de les dades del motor a analitzar. Totes elles expressades
% segons el sistema internacional d'unitats.
M_0=0.85; z=11000; T_SL=288.15; P_SL=101325;
pi_d=0.98; eta_f=0.89; eta_lpc=0.88; eta_hpc=0.86;
pi_b=0.99; eta_b=0.99; eta_hpt=0.91; eta_lpt=0.92;
eta_mH=0.993; eta_mL=0.997; pi_np=0.99; pi_ns=0.99;
T_t4=1450; gamma_c=1.4; gamma_t=1.3; R_g=287; h=43000000; g=9.81;

%% CÀLCULS PREVIS

% Es calculen els coeficients de pressió per l'aire del compressor (Cp_c) i
% de la turbina (Cp_t), segons les gammes depenent de si ha passat o no per
% la cambra de combustió.
Cp_t=(gamma_t/(gamma_t-1))*R_g;
Cp_c=(gamma_c/(gamma_c-1))*R_g;

% S'obtenen també els valors de pressió i temperatura del flux lliure a
% l'altura de funcionament. Amb aquests es calula la velocitat d'entrada de
% l'aire al motor.
[T_0,P_0]=atmosferaISA(z,T_SL,P_SL);
u_0=M_0*sqrt(gamma_c*R_g*T_0);

%% Anàlisi de la relació de compressió del fan

%El set de paràmetres coincideix amb el de la tasca 3.
alpha=10;
pi_f_nominal=1.5737;
pi_lpc=5.3;
pi_hpc=8.5;
OPR=pi_lpc*pi_hpc;

% Es defineix la sensibilitat (variació +- del paràmetre) i els extrems de
% l'interval.
sensibilitat_pi_f=0.5;
pi_f_inf=pi_f_nominal-sensibilitat_pi_f;
pi_f_sup=pi_f_nominal+sensibilitat_pi_f;

% Nombre d'elements del paràmetre variat i definició del vector de relació
% de compressió del fan
N=1000;
vec_pi_f=linspace(pi_f_inf,pi_f_sup,N);

%Definició dels vectors emprats on es recolliran les dades extretes
eta_p=zeros(N,1);
u_9=zeros(N,1);
u_19=zeros(N,1);
M_9_mat=zeros(N,1);
T_9=zeros(N,1);
relacio_area_tovera_primaria=zeros(N,1);
M_19_mat=zeros(N,1);
T_19=zeros(N,1);
relacio_area_tovera_secundaria=zeros(N,1);


%% Estudi de la influència de la relació de compressió del fan

%La funció del for és recollir els resultats desitjats per a tots els valors de
%la relació de compressió del fan anteriorment definits.
  
    for i=1:N
        pi_f=vec_pi_f(i);
        
        % (0)-(2): INLET
        theta_0=1+0.5*(gamma_c-1)*M_0^2; % theta_0=T_t0/T_0
        delta_0=theta_0^(gamma_c/(gamma_c-1)); % delta_0=P_t0/P_0

        % Es calcula la temperatura d'estancament del flux lliure, que roman
        % constant durant la presa d'aire.
        
        T_t0=T_0*theta_0;
        T_t2=T_t0;
        
        % Es calculen també les pressions d'estancament. Aquesta és menor a la presa
        % d'aire degut a l'eficiència d'aquesta.

        P_t0=P_0*delta_0;
        P_t2=P_t0*pi_d;

        % (2)-(2.5): LPC
        % Es calcula la pressió al compressor de baixa pressió gràcies al rendiment
        % isentròpic d'aquest. A partir d'aquest, i de la gamma de l'aire s'obté
        % l'eficiència isentròpica que permet conèixer la temperatura.
        
        P_t25=P_t2*pi_lpc;
        tau_lpc=1+(pi_lpc^((gamma_c-1)/gamma_c)-1)/eta_lpc;
        T_t25=T_t2*tau_lpc;

        % (2.5)-(3): HPC
        % El procediment de càlcul al compressor d'alta pressió és el mateix que al
        % de baixa, partint de la temperatura i la pressió d'aquest últim.
        
        P_t3=P_t25*pi_hpc;
        tau_hpc=1+(pi_hpc^((gamma_c-1)/gamma_c)-1)/eta_hpc;
        T_t3=T_t25*tau_hpc;

        % (3)-(4): COMB
        % Amb les propietats físiques de l'aire a l'entrada i la sortida d'aquesta i
        % l'eficiència isentròpica es pot conèixer la relació de cabal de
        % combustible respecte el cabal d'aire. I amb la relació de pressions
        % d'estancament, la pressió a la sortida.
        f=(Cp_t*T_t4-Cp_c*T_t3)/(eta_b*h-Cp_t*T_t4);
        P_t4=P_t3*pi_b;

        % (4)-(4.5): HPT
        % Per obtenir la relació de temperatures a la turbina d'alta pressió s'usen
        % l'eficiència mecànica en aquesta, així com les temperatures al compressor
        %i a l'entrada de la HPT, tenint en compte la presència de combustible. A
        % partir d'aquesta relació, ja es poden calcular la relació de pressions i
        % la temperatura i pressió a la sortida d'aquest element.

        tau_hpt=1-(1/eta_mH)*(Cp_c/Cp_t)*(1/(1+f))*(T_t25/T_t4)*(tau_hpc-1);
        pi_hpt=(1-(1-tau_hpt)/eta_hpt)^(gamma_t/(gamma_t-1));
        T_t45=T_t4*tau_hpt;
        P_t45=P_t4*pi_hpt;

        % (4.5)-(5): LPT
        %La turbina de baixa pressió segueix el mateix concepte i procediment que
        %la turbina d'alta pressió.
        tau_f=1+(pi_f^((gamma_c-1)/gamma_c)-1)/eta_f;
        tau_lpt=1-(1/eta_mL)*(Cp_c/Cp_t)*(1/(1+f))*(T_t0/T_t4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_f-1));
        pi_lpt=(1-(1-tau_lpt)/eta_lpt)^(gamma_t/(gamma_t-1));
        T_t5=T_t45*tau_lpt;
        P_t5=P_t45*pi_lpt;

        % (5)-(9): PRIMARY NOZZLE
        %Les temperatures d'estancament de la tovera primària es mantenen
        %constants, mentre que la pressió disminueix lleugerament a causa de
        %l'eficiència de la tovera.
        
        P_t9=P_t5*pi_np;
        T_t9=T_t5;
        
        %Relació de pressió d'estancament a la tovera i a l'entrada del motor.
        pi_primari_total=P_t9/P_0;

        %En aquest cas si cal imposar tovera adaptada degut a que amb la
        %variació de la relació de compressió del fan si que hi ha moments
        %on la la tovera primària està ofegada.
        
        P_9=P_0;
        pi_primari=1;
        M_9_mat(i)=sqrt(2/(gamma_t-1)*(pi_primari_total^((gamma_t-1)/gamma_t)-1));
        
        T_9(i)=T_t9/(1+0.5*(gamma_t-1)*M_9_mat(i)^2);
        u_9(i)=M_9_mat(i)*sqrt(gamma_t*R_g*T_9(i));

        M_8_mat=1;
        MFP_8=sqrt(gamma_t)*(M_8_mat/(1+0.5*(gamma_t-1)*M_8_mat^2)^((gamma_t+1)/(2*(gamma_t-1))));
        MFP_9=sqrt(gamma_t)*(M_9_mat(i)/(1+0.5*(gamma_t-1)*M_9_mat(i)^2)^((gamma_t+1)/(2*(gamma_t-1))));
        relacio_area_tovera_primaria(i)=MFP_8/MFP_9;
        
        % (0)-(1.3): SECONDARY INLET
        % A la presa secondària d'aire, la temperatura d'estancament és major que la
        % del flux lliure.
        T_t13=T_0*theta_0*tau_f;

        % (1.3)-(1.9): SECONDARY NOZZLE
        % A la presa secondària d'aire, la temperatura d'estancament és major que la
        % del flux lliure.
        T_t19=T_t13;
        
        % pi_secundari_total=P_t19/P_0
        pi_secundari_total=delta_0*pi_d*pi_f*pi_ns; % pi_secundari_total=P_t19/P_0
        
        %De la mateixa manera que a la tasca 4, s'imposa tovera adaptada.
        M_19_mat(i)=sqrt(2/(gamma_c-1)*(pi_secundari_total^((gamma_c-1)/gamma_c)-1));
        P_19=P_0;
        pi_secundari=1;
        
        T_19(i)=T_t19/(1+(gamma_c-1)/2*M_19_mat(i)^2);
        u_19(i)=M_19_mat(i)*sqrt(gamma_c*R_g*T_19(i));
        
        %Càlcul dels mass flow parameters i de la relació d'àrees 
        M_18_mat=1;
        MFP_18=sqrt(gamma_t)*(M_18_mat/(1+0.5*(gamma_c-1)*M_18_mat^2)^((gamma_c+1)/(2*(gamma_c-1))));
        MFP_19=sqrt(gamma_t)*(M_19_mat(i)/(1+0.5*(gamma_c-1)*M_19_mat(i)^2)^((gamma_c+1)/(2*(gamma_c-1))));
        relacio_area_tovera_secundaria(i)=MFP_18/MFP_19;

        % Càlcul del Thrust específic
        Specific_Thrust_primari=(1+f)*u_9(i)-u_0+(1+f)*T_9(i)/u_9(i)*R_g*(1-1/pi_primari);
        Specific_Thrust_secundari=alpha*(u_19(i)-u_0)+alpha*(T_19(i)*R_g)/u_19(i)*(1-1/pi_secundari);
        Specific_Thrust=Specific_Thrust_primari+Specific_Thrust_secundari;

        % Specific impuls
        I_sp=Specific_Thrust/(f*g);

        eta_p(i)=2*Specific_Thrust*u_0/((1+f)*u_9(i)^2+alpha*u_19(i)^2-(1+alpha)*u_0^2);   
        
        if imag(eta_p(i))~=0
            eta_p(i)=0;
        end
        if imag(u_9(i))~=0
            u_9(i)=0;
        end
        if imag(u_19(i))~=0
            u_19(i)=0;
        end
        if imag(M_9_mat(i))~=0
            M_9_mat(i)=0;
        end
        if imag(relacio_area_tovera_primaria(i))~=0
            relacio_area_tovera_primaria(i)=0;
        end
        if imag(M_19_mat(i))~=0
            M_19_mat(i)=0;
        end
        if imag(relacio_area_tovera_secundaria(i))~=0
            relacio_area_tovera_secundaria(i)=0;
        end
        
    end
    
    %% Gràfiques
    
    figure(1);
    set(figure(1),'Renderer', 'painters', 'Position', [20 490 550 450]);
    plot(vec_pi_f,u_9)
    ylabel('Velocitat de sortida (flux primari) $u_9\;\mathrm{(m/s)}$','Interpreter','latex');
    xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
    grid on
    grid minor
    ax = gca;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha=0.3;
    ax.MinorGridColor = [0, 0, 0];
    ax.MinorGridAlpha=0.3;
    xlim([vec_pi_f(1),vec_pi_f(end)]);
    
    figure(2);
    set(figure(2),'Renderer', 'painters', 'Position', [20+600 490 550 450]);
    plot(vec_pi_f,u_19)
    ylabel('Velocitat de sortida (flux secundari) $u_{19}\;\mathrm{(m/s)}$','Interpreter','latex');
    xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
    grid on
    grid minor
    ax = gca;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha=0.3;
    ax.MinorGridColor = [0, 0, 0];
    ax.MinorGridAlpha=0.3;
    xlim([vec_pi_f(1),vec_pi_f(end)]);
    
    figure(3);
    set(figure(3),'Renderer', 'painters', 'Position', [20+1200 490 550 450]);
    plot(vec_pi_f,eta_p,'r')
    ylabel('Rendiment propulsiu $\eta_p$','Interpreter','latex');
    xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
    grid on
    grid minor
    ax = gca;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha=0.3;
    ax.MinorGridColor = [0, 0, 0];
    ax.MinorGridAlpha=0.3;
    xlim([vec_pi_f(1),vec_pi_f(end)]);
    
% Gràfiques addicionals per l'annex

%     figure(4);
%     set(figure(4),'Renderer', 'painters', 'Position', [20 490 400 350]);
%     yyaxis left;
%     plot(vec_pi_f,M_9_mat)
%     ylabel('N\''umero de Mach  $M_{9,\,mat}$ (tovera adaptada)','Interpreter','latex');
%     xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
%     yyaxis right;
%     plot(vec_pi_f,relacio_area_tovera_primaria)
%     ylabel('Relaci\''o d''\`arees $A_{9}/A_{8}$ (tovera adaptada)','Interpreter','latex');
%     legend('Mach $M_{9,\,mat}$','Relaci\''o $A_{9}/A_{8}$','Location','northwest','Interpreter','latex');
%     grid on
%     grid minor
%     ax = gca;
%     ax.GridColor = [0, 0, 0];
%     ax.GridAlpha=0.3;
%     ax.MinorGridColor = [0, 0, 0];
%     ax.MinorGridAlpha=0.3;
%     xlim([vec_pi_f(1),vec_pi_f(end)]);
%     
%         figure(5);
%     set(figure(5),'Renderer', 'painters', 'Position', [20+400 490 400 350]);
%     plot(vec_pi_f,T_9,'r')
%     ylabel('Temperatura $T_{9}$ ($\mathrm{K}$)','Interpreter','latex');
%     xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
%     grid on
%     grid minor
%     ax = gca;
%     ax.GridColor = [0, 0, 0];
%     ax.GridAlpha=0.3;
%     ax.MinorGridColor = [0, 0, 0];
%     ax.MinorGridAlpha=0.3;
%     xlim([vec_pi_f(1),vec_pi_f(end)]);
%     
%         figure(6);
%     set(figure(6),'Renderer', 'painters', 'Position', [20+800 490 400 350]);
%     yyaxis left;
%     plot(vec_pi_f,M_19_mat)
%     ylabel('N\''umero de Mach  $M_{19,\,mat}$ (tovera adaptada)','Interpreter','latex');
%     xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
%     yyaxis right;
%     plot(vec_pi_f,relacio_area_tovera_secundaria)
%     ylabel('Relaci\''o d''\`arees $A_{19}/A_{18}$ (tovera adaptada)','Interpreter','latex');
%     xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
%     legend('Mach $M_{19,\,mat}$','Relaci\''o $A_{19}/A_{18}$','Location','northwest','Interpreter','latex')
%     grid on
%     grid minor
%     ax = gca;
%     ax.GridColor = [0, 0, 0];
%     ax.GridAlpha=0.3;
%     ax.MinorGridColor = [0, 0, 0];
%     ax.MinorGridAlpha=0.3;
%     xlim([vec_pi_f(1),vec_pi_f(end)]);
%     
%             figure(7);
%     set(figure(7),'Renderer', 'painters', 'Position', [20+1200 490 400 350]);
%     plot(vec_pi_f,T_19,'r')
%     ylabel('Temperatura $T_{19}$ ($\mathrm{K}$)','Interpreter','latex');
%     xlabel('Ratio de pressions del fan $\pi_{f}$','Interpreter','latex')
%     grid on
%     grid minor
%     ax = gca;
%     ax.GridColor = [0, 0, 0];
%     ax.GridAlpha=0.3;
%     ax.MinorGridColor = [0, 0, 0];
%     ax.MinorGridAlpha=0.3;
%     xlim([vec_pi_f(1),vec_pi_f(end)]);





