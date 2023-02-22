clear

%% Drop-Wall Benchmark (Movie S7 From de Ruiter et al. 2014)
R=0.69e-3;
U=0.26;
surf_tens=19.7e-3;
g=9.81;
rho_l=913.0;
mu_l=5.0e-3;
mu_v=1.827e-5;
A_H_dim=3.7e-20;
mfp=69e-9; % Mean free path at atmospheric pressure
press_r=1.0; % ratio to atmospheric pressure
use_gke=true;
drop_drop=false;
folder='RESLT_benchmark_drop_wall';
filename='benchmark_drop_wall_run.sh';
[Re,Grav,Cap,VisR,Ham,Kn_R]=compute_parameters(R,U,rho_l,mu_l,g,surf_tens,mu_v,A_H_dim,mfp,press_r);
output_script(filename,folder,Re,Grav,Cap,VisR,Ham,Kn_R,use_gke,drop_drop);

%% Drop-Drop Benchmark (Fig 3 from Pan et al 2008)
R=167.6e-6;
U=0.496;
rho_l =762.0;
mu_l=2.128e-3;
g=0.0;
surf_tens=2.65e-2;
mu_v=1.827e-5;
A_H_dim=5.0e-20;
mfp=69e-9; % Mean free path at atmospheric pressure
press_r=1.0; % ratio to atmospheric pressure
use_gke=true;
drop_drop=true;
folder='RESLT_benchmark_drop_drop';
filename='benchmark_drop_drop_run.sh';
[Re,Grav,Cap,VisR,Ham,Kn_R]=compute_parameters(R,U,rho_l,mu_l,g,surf_tens,mu_v,A_H_dim,mfp,press_r);
output_script(filename,folder,Re,Grav,Cap,VisR,Ham,Kn_R,use_gke,drop_drop);

%% Tetradecane Example
R=150.0e-6;
rho_l=762.0;
mu_l=2.128e-3;
g=0;
surf_tens=2.65e-2;
mu_v=1.84e-5;
A_H_dim=5.0e-20;
mfp=69e-9; % Mean free path at atmospheric pressure
press_r=1.0; % ratio to atmospheric pressure
use_gke=true;
drop_drop=true;
We=12.0;
U=sqrt(We*surf_tens/(4.0*R*rho_l));

folder='RESLT_example';
filename='example_run.sh';
[Re,Grav,Cap,VisR,Ham,Kn_R]=compute_parameters(R,U,rho_l,mu_l,g,surf_tens,mu_v,A_H_dim,mfp,press_r);
output_script(filename,folder,Re,Grav,Cap,VisR,Ham,Kn_R,use_gke,drop_drop);


