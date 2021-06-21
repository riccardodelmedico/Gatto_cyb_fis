%% Dati su vaccinazione in Italia (al 2021/06/9)
T=readtable('somministrazioni-vaccini-latest.csv');

pop= 59257566; %istat 2020
n=height(T);
vaccini=[T(1:n,1), T(1:n,7), T(1:n,8)];
vaccini_matrix_timetable= table2timetable(vaccini);
vaccini_daily= retime(vaccini_matrix_timetable, 'daily', 'sum');
vaccini_table= timetable2table(vaccini_daily);
vaccini_matrix= table2array(vaccini_table(:,2:3));
prima_dose= vaccini_matrix(:,1);
seconda_dose=vaccini_matrix(:,2);

prima_dose_norm= prima_dose./pop;
seconda_dose_norm= seconda_dose./pop;

clear n vaccini vaccini_daily vaccini_matrix vaccini_matrix vaccini_matrix_timetable ...
    vaccini_table prima_dose seconda_dose T