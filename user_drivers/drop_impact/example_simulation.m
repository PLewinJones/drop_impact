clear

%% Tetradecane
R=150.0e-6;
U=0.5;
rho_l=762.0;
mu_l=2.128e-3;

surf_tens=2.65e-2;
mu_v=1.84e-5;
A_H_dim=5.0e-20;
mfp=69e-9;
g=0;
press_r=1.0;
use_gke=true;

T_end=10.0; %The simulation is cut off if the drop bounces, so this just needs to be large enough.

folder='RESLT_example';
Re=U*R*rho_l/mu_l;
Grav=g*R^2*rho_l/(U*mu_l);
Cap=mu_l*U/surf_tens;
VisR=mu_v/mu_l;
Ham=A_H_dim/(6*pi*mu_l*U*R^2);
Lam=mfp/(R*press_r);

options=append('--reynolds ',num2str(Re,8),...
  ' --gravity ',num2str(Grav,8),...
  ' --capilliary ',num2str(Cap,8),...
  ' --viscosity ',num2str(VisR,8),...
  ' --hamaker ',num2str(Ham,8),...
   ' --tmax ',num2str(T_end,8));
if(use_gke)
  options=append(options,' --mfp ',num2str(Lam,8),' --usegke ');
end
options=append(options,' --dropdrop ');
filename=append('run_drop_impact.sh');
fileid=fopen(filename,'w');
fprintf(fileid,'#!/bin/bash \n');
fprintf(fileid,'\n');
fprintf(fileid,append('folder=',folder,' \n'));
fprintf(fileid,'logfile=$folder"/oomph-out.txt"\n');
fprintf(fileid,append('parameters="',options,'"\n'));
fprintf(fileid,'if [[ -d $folder ]]; then\n');
fprintf(fileid,'echo "$folder already exists, stopping simulation"\n');
fprintf(fileid,'exit\n');
fprintf(fileid,'fi\n');
fprintf(fileid,'mkdir $folder\n');
fprintf(fileid,'make\n');
fprintf(fileid,'./drop_impact --folder $folder $parameters &> $logfile\n');
fclose(fileid);


