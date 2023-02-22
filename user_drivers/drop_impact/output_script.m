function output_script(filename,folder,Re,Grav,Cap,VisR,Ham,Kn_R,use_gke,drop_drop)

T_end=10.0; %The simulation is cut off if the drop bounces, so this just needs to be large enough.

options=append('--reynolds ',num2str(Re,8),...
  ' --gravity ',num2str(Grav,8),...
  ' --capilliary ',num2str(Cap,8),...
  ' --viscosity ',num2str(VisR,8),...
  ' --hamaker ',num2str(Ham,8),...
   ' --tmax ',num2str(T_end,8));
if(use_gke)
options=append(options,' --mfp ',num2str(Kn_R,8),' --usegke ');
end
if(drop_drop)
options=append(options,' --dropdrop ');
end

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

end