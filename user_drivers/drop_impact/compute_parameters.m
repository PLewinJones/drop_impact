function [Re,Grav,Cap,VisR,Ham,Kn_R]=...
  compute_parameters(R,U,rho_l,mu_l,g,surf_tens,mu_v,A_H_dim,mfp,press_r)
Re=U*R*rho_l/mu_l;
Grav=g*R^2*rho_l/(U*mu_l);
Cap=mu_l*U/surf_tens;
VisR=mu_v/mu_l;
Ham=A_H_dim/(6*pi*mu_l*U*R^2);
Kn_R=mfp/(R*press_r);
end