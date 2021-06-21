%% Ricerca del tempo di raddoppio con R0=3.6
global lambda x0 N R0 Lvect 

format longg
R0=3.6;
parameters_vaccini;
dati_vaccini;
E0= 10/pop;
Lvect = zeros(N,1);
x0= [1-1*E0; 1*E0; zeros(22,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, 30, x0); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30, 24);
for j=1:1:30
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,5)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1)) %calcoli nel report, è il coefficiente angolare
t_radd = log(2)/tasso_expo %tempo di raddoppio


%% ritaratura R0, mettendo il tempo di raddoppio a CIRCA 3, A ''OCCHIO''
R0 = 2.65;
r0_raddoppio;

Lvect = zeros(N,1);
x0= [1-1*E0; 1*E0; zeros(22,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, 30, x0); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30, 24);
for j=1:1:30
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,5)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% Ricerca del tempo di raddoppio con R0=3.6 (waterfall ODE su E e H)
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) rifacciamo lo studio
R0=3.6;
r0_raddoppio;

Lvect = zeros(N,1);
x0_casc = [1-1*E0; 1*E0; zeros(29,1)]; 

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascate, 0, 1, 30, x0_casc); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,31);
for j=1:1:30
    for i= 1:1:31
        soluzione(j,i)=x_vaccini_tot(31*(j-1)+i);
    end
end

E = [soluzione(:,3)];
P = [soluzione(:,4)];
I = [soluzione(:,5)];
A = [soluzione(:,6)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% Ritaratura di R0 (waterfall ODE E,H)
R0=2.65;
r0_raddoppio;

Lvect = zeros(N,1);
x0_casc = [1-1*E0; 1*E0; zeros(29,1)]; 

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascate, 0, 1, 30, x0_casc); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,31 );
for j=1:1:30
    for i= 1:1:31
        soluzione(j,i)=x_vaccini_tot(31*(j-1)+i);
    end
end

E = [soluzione(:,3)];
P = [soluzione(:,4)];
I = [soluzione(:,5)];
A = [soluzione(:,6)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% Ricerca del tempo di raddoppio con R0=3.6 (waterfall ODE su E, I, A)
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 
R0=3.6;
r0_raddoppio;

Lvect = zeros(N,1);
x0_casc_inf = [1-1*E0; 1*E0; zeros(40,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, 30, x0_casc_inf); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,42);
for j=1:1:30
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,7)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% Ritaratura di R0
R0 = 2.65;
r0_raddoppio;

Lvect = zeros(N,1);
x0_casc_inf = [1-1*E0; 1*E0; zeros(40,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, 30, x0_casc_inf); 
toc

soluzione= zeros(30,42);
for j=1:1:30
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,7)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo
