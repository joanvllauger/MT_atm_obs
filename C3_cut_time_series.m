addpath './aux_functions/'
% Cridam C1_diagnostic per diferents esdeveniments
%--------------------------------------------------------------------------
close all
clear all
T=readtable('List_of_events.xlsx');
times=datenum(table2array(T(:,1)));
t_in=datenum(table2array(T(:,2)));
t_out=datenum(table2array(T(:,3)));

path_mat='./data/';

%%
%-- path dades
arx_dades='Atm_pres_all.mat';
disp('Pressió atmosferica')

%--- Variables interesants que depenen del temps
variables={'mareografs','lon','lat','time','Patm'};

for nt=1:length(times)
    disp('--------------------------------------------------------')
    disp(datestr(times(nt)))
    t_ref=times(nt);
    tini=t_in(nt)-3/24;
    tfi=t_out(nt)+3/24;

    % 0- Importam i netejam les dades
    [dades,b]=F1_clean_ST_v2(path_mat,arx_dades,variables,tini,tfi,0);

    for nv=1:length(variables)
        eval([variables{nv},'=dades.',variables{nv},';'])
    end

    disp('Dades importades')

    % Cream una carpeta per guardar les dades
    folder=['./meteotsunamis/meteotsunami_',datestr(t_ref,'yyyymmdd')];
    mkdir(folder)


    %--------------------------------------------------------------------------
    save([folder,'/atm_info.mat'],variables{:},'-v7.3');
    disp('Guardat')
end
