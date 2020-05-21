% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

%En aquest codi, es realitzen les gràfiques de les prestacions del turbofan
%en funció de cadascun del paràmetres no definits, al voltant del valor
%nominal de l'apartat 1 i dintre d'un rang realista.

clear
close all
clc


% Activar la següent variable com a 'true' si es desitja visualitzar les
% contribucions de cada flux, primari i secundari, a les gràfiques de
% l'empenyiment específic segons la relació de compressió del fan
mostrar_contribucions=false;

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

%% DEFINICIÓ DE RANGS D'ESTUDI

% El grup de paràmetres escollits a la tasca 1 s'estableix com a nominal
alpha_nominal=8;
pi_f_nominal=1.6;
pi_lpc_nominal=4.3;
pi_hpc_nominal=7.5;
valors_nominals=[alpha_nominal,pi_f_nominal,pi_lpc_nominal,pi_hpc_nominal];

% Es defineix la sensibilitat (variació +- del paràmetre)
sensibilitat_alpha=2;
sensibilitat_pi_f=0.5;
sensibilitat_pi_lpc=1;
sensibilitat_pi_hpc=1;

% Nombre d'elements per paràmetre variat
N=1000;
valors_Specific_Thrust=zeros(N,1);
valors_I_sp=zeros(N,1);
Specific_Thrust_primari=zeros(N,1);
Specific_Thrust_secundari=zeros(N,1);

% Definició de límits per a la relació de derivació
alpha_inf=alpha_nominal-sensibilitat_alpha;
alpha_sup=alpha_nominal+sensibilitat_alpha;
% Definició del vector de valors de la relació de derivació
vec_alpha=linspace(alpha_inf,alpha_sup,N);

% Definició de límits per a la relació de pressions del fan
pi_f_inf=pi_f_nominal-sensibilitat_pi_f;
pi_f_sup=pi_f_nominal+sensibilitat_pi_f;
% Definició del vector de valors de la relació de pressions del fan
vec_pi_f=linspace(pi_f_inf,pi_f_sup,N);

% Definició de límits per a la relació de pressions del LPC
pi_lpc_inf=pi_lpc_nominal-sensibilitat_pi_lpc;
pi_lpc_sup=pi_lpc_nominal+sensibilitat_pi_lpc;
% Definició del vector de valors de la relació de pressions del LPC
vec_pi_lpc=linspace(pi_lpc_inf,pi_lpc_sup,N);

% Definició de límits per a la relació de pressions del HPC
pi_hpc_inf=pi_hpc_nominal-sensibilitat_pi_hpc;
pi_hpc_sup=pi_hpc_nominal+sensibilitat_pi_hpc;
% Definició del vector de valors de la relació de pressions del HPC
vec_pi_hpc=linspace(pi_hpc_inf,pi_hpc_sup,N);

% Es composa una matriu amb tots els vectors de valors
matriu_vectors_variables=[vec_alpha;vec_pi_f;vec_pi_lpc;vec_pi_hpc];


%% Estudi de la influència de cada paràmetre

% Es crea un bucle que analitza de forma individual la influència de cada
% paràmetre
for j=1:4
    
    % En funció de la variable, el format de les gràfiques serà diferent
    % A més, es fixen els paràmetres que no canvien
    if j==1     % Cas de la relació de derivació
       pi_f=pi_f_nominal;
       pi_lpc=pi_lpc_nominal;
       pi_hpc=pi_hpc_nominal;
       noms_variable='\textit{By-pass ratio} $\alpha$';
    elseif j==2 % Cas de la ràtio de pressions del fan
       alpha=alpha_nominal;
       pi_lpc=pi_lpc_nominal;
       pi_hpc=pi_hpc_nominal;
       noms_variable='R\`atio de pressions del fan $\pi_{f}$';
    elseif j==3 % Cas de la ràtio de pressions del LPC
       alpha=alpha_nominal;
       pi_f=pi_f_nominal;
       pi_hpc=pi_hpc_nominal;
       noms_variable='R\`atio de pressions del LPC $\pi_{LPC}$';
    else % Cas de la ràtio de pressions del HPC
       alpha=alpha_nominal;
       pi_f=pi_f_nominal;
       pi_lpc=pi_lpc_nominal;
       noms_variable='R\`atio de pressions del HPC $\pi_{HPC}$';
    end 
    
    % S'inicia un bucle anidat que recorre tot el rang de valors del
    % paràmetre variable
    for i=1:N 
        % S'obté el vector de valors del paràmetre variable
        if j==1
            alpha=matriu_vectors_variables(j,i);
        elseif j==2
            pi_f=matriu_vectors_variables(j,i);
        elseif j==3
            pi_lpc=matriu_vectors_variables(j,i);
        else
            pi_hpc=matriu_vectors_variables(j,i);
        end
        
        % Resolució del motor d'acord a la metodologia analítica de la
        % tasca 1 per a cada valor del paràmetre variable
        
        % (0)-(2): INLET
        theta_0=1+0.5*(gamma_c-1)*M_0^2; % theta_0=T_t0/T_0
        delta_0=theta_0^(gamma_c/(gamma_c-1)); % delta_0=P_t0/P_0

        T_t0=T_0*theta_0;
        T_t2=T_t0;

        P_t0=P_0*delta_0;
        P_t2=P_t0*pi_d;

        % (2)-(2.5): LPC
        P_t25=P_t2*pi_lpc;
        tau_lpc=1+(pi_lpc^((gamma_c-1)/gamma_c)-1)/eta_lpc;
        T_t25=T_t2*tau_lpc;

        % (2.5)-(3): HPC
        P_t3=P_t25*pi_hpc;
        tau_hpc=1+(pi_hpc^((gamma_c-1)/gamma_c)-1)/eta_hpc;
        T_t3=T_t25*tau_hpc;

        % (3)-(4): COMB
        f=(Cp_t*T_t4-Cp_c*T_t3)/(eta_b*h-Cp_t*T_t4);
        P_t4=P_t3*pi_b;

        % (4)-(4.5): HPT
        tau_hpt=1-(1/eta_mH)*(Cp_c/Cp_t)*(1/(1+f))*(T_t25/T_t4)*(tau_hpc-1);
        pi_hpt=(1-(1-tau_hpt)/eta_hpt)^(gamma_t/(gamma_t-1));
        T_t45=T_t4*tau_hpt;
        P_t45=P_t4*pi_hpt;

        % (4.5)-(5): LPT
        tau_f=1+(pi_f^((gamma_c-1)/gamma_c)-1)/eta_f;
        tau_lpt=1-(1/eta_mL)*(Cp_c/Cp_t)*(1/(1+f))*(T_t0/T_t4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_f-1));
        pi_lpt=(1-(1-tau_lpt)/eta_lpt)^(gamma_t/(gamma_t-1));
        T_t5=T_t45*tau_lpt;
        P_t5=P_t45*pi_lpt;

        % (5)-(9): PRIMARY NOZZLE
        P_t9=P_t5*pi_np;
        T_t9=T_t5;
        pi_primari_total=P_t9/P_0; %Poner esto o cadena?
        M_9=sqrt(2/(gamma_t-1)*((pi_primari_total)^((gamma_t-1)/gamma_t)-1));
        if M_9<1
            P_9=P_0;
            pi_primari=1;
        else
            M_9=1;
            pi_primari=pi_primari_total*(2/(gamma_t+1))^(gamma_t/(gamma_t-1)); % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que és fluxe secundari
            P_9=pi_primari*P_0;
        end
        T_9=T_t9/(1+0.5*(gamma_t-1)*M_9^2);
        u_9=M_9*sqrt(gamma_t*R_g*T_9);

        % (0)-(1.3): SECONDARY INLET
        T_t13=T_0*theta_0*tau_f;

        % (1.3)-(1.9): SECONDARY NOZZLE
        T_t19=T_t13;

        pi_secundari_total=delta_0*pi_d*pi_f*pi_ns; % pi_secundari_total=P_t19/P_0
        M_19=sqrt(2/(gamma_c-1)*((pi_secundari_total)^((gamma_c-1)/gamma_c)-1));
        if M_19<1
            P_19=P_0;
            pi_secundari=1;
            P_t19=pi_secundari_total*P_0;
        else
            M_19=1;
            pi_secundari=pi_secundari_total*(2/(gamma_c+1))^(gamma_c/(gamma_c-1)); % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que és fluxe secundari
            P_19=pi_secundari*P_0;
            P_t19=pi_secundari_total*P_0;
        end
        T_19=T_t19/(1+(gamma_c-1)/2*M_19^2);
        u_19=M_19*sqrt(gamma_c*R_g*T_19);

        % Càlcul de l'empenta específica
        Specific_Thrust_primari(i)=(1+f)*u_9-u_0+(1+f)*T_9/u_9*R_g*(1-1/pi_primari);
        Specific_Thrust_secundari(i)=alpha*(u_19-u_0)+alpha*(T_19*R_g)/u_19*(1-1/pi_secundari);
        Specific_Thrust=Specific_Thrust_primari(i)+Specific_Thrust_secundari(i);  
        
        % Càlcul de l'impuls específic
        I_sp=Specific_Thrust/(f*g);
        
        % Degut a la possibilitat d'obtenir valors complexos en el resultat
        % a causa de la impossibilitat de complir el balanç de potències,
        % els resultats es consideran nuls en aquesta situació (són físicament impossibles)
        if imag(I_sp)~=0||imag(Specific_Thrust)~=0
        Specific_Thrust=0;
        I_sp=0;
        u_9=0;
        u_19=0;
        f=0;
        end

        % Emmagatzemament de valors associats a la variació d'un dels
        % paràmetres
        valors_Specific_Thrust(i)=Specific_Thrust;
        valors_I_sp(i)=I_sp;
        
        
              
    end
    
    % Representacions gràfiques associades a la variació d'aquest paràmetre
    
    % Gràfica de l'impuls específic
    figure(2*j-1);
    set(figure(2*j-1),'Renderer', 'painters', 'Position', [20+(j-1)*390 490 385 300]);
    plot(matriu_vectors_variables(j,:),valors_I_sp,'b')
    ylabel('Impuls espec\''ific $I_{sp}\;\mathrm{(s)}$','Interpreter','latex');
    xlabel(noms_variable,'Interpreter','latex')
    grid on
    grid minor
    ax = gca;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha=0.3;
    ax.MinorGridColor = [0, 0, 0];
    ax.MinorGridAlpha=0.3;
    xlim([matriu_vectors_variables(j,1),matriu_vectors_variables(j,end)]);
    if(j==2)
    xticks([pi_f_inf:0.1:pi_f_sup])
    elseif(j==3)
    xticks([pi_lpc_inf:0.5:pi_lpc_sup])
    end
    
    % Gràfica de l'empenta específica
    figure(2*j);
    set(figure(2*j),'Renderer', 'painters', 'Position', [20+(j-1)*390 100 385 300]);
    plot(matriu_vectors_variables(j,:),valors_Specific_Thrust,'r');
    if mostrar_contribucions==true
    hold on
    plot(matriu_vectors_variables(j,:),Specific_Thrust_primari,'m--',matriu_vectors_variables(j,:),Specific_Thrust_secundari,'r--');
    legend('Contribuci\''o total dels fluxos a $F/\dot{m}$','Contribuci\''o del flux primari','Contribuci\''o del flux secundari','Location','southwest','Interpreter','latex')
    end
    ylabel('Empenyiment espec\''ific $F/\dot{m}\;\mathrm{(m/s)}$','Interpreter','latex')
    xlabel(noms_variable,'Interpreter','latex')
    grid on
    grid minor
    ax = gca;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha=0.3;
    ax.MinorGridColor = [0, 0, 0];
    ax.MinorGridAlpha=0.3;
    xlim([matriu_vectors_variables(j,1),matriu_vectors_variables(j,end)]);
    if(j==2)
    xticks([pi_f_inf:0.1:pi_f_sup])
    elseif(j==3)
    xticks([pi_lpc_inf:0.5:pi_lpc_sup])


    end
   
end
