// Interface lubrication element. Unified problems
#include "surface_lubrication_elements.h"


namespace oomph
{
  // Define the default physical value to be one
  double InterfaceLubricationElement::Default_Physical_Constant_Value = 1.0;

  //=====================================================================
  /// Get the Vapour pressure
  //=====================================================================
  double InterfaceLubricationElement::interpolated_Pv(const Vector<double>& s)
  {
    // Find number of nodes
    unsigned n_node = this->nnode();

    // Local shape function
    Shape psi(n_node);

    // Find values of shape function
    this->shape(s, psi);

    // Initialise value of Pv
    double Pv = 0.0;

    // Loop over the local nodes and sum
    for (unsigned l = 0; l < n_node; l++)
    {
      Pv += this->nodal_value(l, this->Pv_index[l]) * psi(l);
    }

    return (Pv);
  }

  //=====================================================================
  /// Get the Vapour pressure
  //=====================================================================
  double InterfaceLubricationElement::interpolated_Al(const Vector<double>& s)
  {
    // Find number of nodes
    unsigned n_node = this->nnode();

    // Local shape function
    Shape psi(n_node);

    // Find values of shape function
    this->shape(s, psi);

    // Initialise value of Pv
    double Al = 0.0;

    // Loop over the local nodes and sum
    for (unsigned l = 0; l < n_node; l++)
    {
      Al += this->nodal_value(l, this->Al_index[l]) * psi(l);
    }

    return (Al);
  }
  //=======================================================================
  /// Overload the Helper function to calculate the residuals and
  /// jacobian entries. This particular function ensures that the
  /// additional entries are calculated inside the integration loop
  void InterfaceLubricationElement::
    add_additional_residual_contributions_interface(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag,
      const Shape& psif,
      const DShape& dpsifds,
      const DShape& dpsifdS,
      const DShape& dpsifdS_div,
      const Vector<double>& s,
      const Vector<double>& interpolated_x,
      const Vector<double>& interpolated_n,
      const double& W,
      const double& J)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Storage for the local equation numbers and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Surface Thin-film equation

    // Find the index at which the velocity is stored
    Vector<unsigned> u_index = this->U_index_interface;

    // Read out the physical constants
    const double Vis_r = this->vapour_viscosity_ratio();
    const double Ha = this->hamaker();
    const double Mfp = this->mean_free_path();
    const double FilmCutoff = this->film_cutoff();

    bool DD_flag = drop_drop_flag();
    double C_DD = 0.0;
    if (DD_flag == true)
    {
      C_DD = 1.0;
    }

    // Get the GKE factors

    double GKE_Tau_P = 1.0;
    double GKE_Tau_C = 1.0;
    double GKE_Phi_P = 1.0;
    double GKE_Phi_C = 1.0;
    double Kn = Mfp /( (1.0+C_DD)*interpolated_x[1]);

    get_GKE_Tau_P(Kn, GKE_Tau_P);
    get_GKE_Tau_C(Kn, GKE_Tau_C);
    get_GKE_Phi_P(Kn, GKE_Phi_P);
    get_GKE_Phi_C(Kn, GKE_Phi_C);

    

    // Now calculate the vapour pressure at this point
    // Assuming the same shape functions are used (which they are)
    double interpolated_Pv = 0.0;

    // The mesh and fluid velocity, and grad of vapour pressure
    const unsigned n_dim = this->node_pt(0)->ndim();
    Vector<double> interpolated_u(n_dim, 0.0);
    Vector<double> interpolated_grad_Pv(n_dim, 0.0);
    Vector<double> interpolated_grad_Al(n_dim, 0.0);
    Vector<double> interpolated_grad_r(n_dim, 0.0);
    Vector<double> interpolated_grad_h(n_dim, 0.0);
    // Vector<double> interpolated_mesh_velocity(n_dim, 0.0);

    double interpolated_dhdr = 0.0;
    double dhdt = 0.0;

    // Are we on the base of the drop?
    bool film_flag = false;
    if (interpolated_n[1] < -FilmCutoff)
    {
      film_flag = true;
    }

    // Loop over the shape functions
    for (unsigned l = 0; l < n_node; l++)
    {
      const double psi = psif(l);
      const double Pv_ = this->nodal_value(l, this->Pv_index[l]);
      const double Al_ = this->nodal_value(l, this->Al_index[l]);
      const double r_ = this->nodal_position(l, 0);
      const double h_ = this->nodal_position(l, 1);
      interpolated_Pv += Pv_ * psi;
      // Velocity and Mesh Velocity
      for (unsigned j = 0; j < n_dim; j++)
      {
        const double u_ = this->nodal_value(l, u_index[j]);
        interpolated_u[j] += u_ * psi;
        interpolated_grad_Pv[j] += Pv_ * dpsifdS(l, j);
        interpolated_grad_Al[j] += Al_ * dpsifdS(l, j);
        interpolated_grad_r[j] += r_ * dpsifdS(l, j);
        interpolated_grad_h[j] += h_ * dpsifdS(l, j);
      }
    }


    if (std::abs(interpolated_grad_r[0]) >
        0.0) // Stop div by 0, should be where p_v is pinned anyway
    {
      interpolated_dhdr += interpolated_grad_h[0] / interpolated_grad_r[0];
    }


    dhdt = interpolated_u[1] - (interpolated_u[0]) * interpolated_dhdr;


    // Now we add the various terms to the appropriate residuals and jacobians
    // The main loop is over the nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      // Read out the apprporiate local equation for the Pv residual and
      // Jacobian terms
      local_eqn = this->nodal_local_eqn(l, this->Pv_index[l]);

      // If not a boundary condition
      if (local_eqn >= 0)
      {
        // Add Residual terms
        // First part
        residuals[local_eqn] +=
          dhdt  * psif(l) * W * J;
        // Second Part
        double velocity_term = 0.0;

        velocity_term +=
          interpolated_x[1] *
          ((interpolated_grad_Pv[0] / interpolated_grad_r[0]) * GKE_Phi_P *
             interpolated_x[1] * interpolated_x[1]*(1.0+3.0*C_DD) / 12.0 -
           GKE_Phi_C * interpolated_u[0]*(1.0+C_DD) / 2.0) *
          dpsifdS(l, 0) / interpolated_grad_r[0];


        residuals[local_eqn] += velocity_term * W * J;

        // We also need to compute the jacobian terms
        if (flag)
        {
          // Loop over the nodes again, this time for the unknowns
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            // Get the unknown pv_index
            local_unknown = this->nodal_local_eqn(l2, this->Pv_index[l2]);
            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                (GKE_Phi_P * interpolated_x[1] * interpolated_x[1] *
                 interpolated_x[1]*(1.0+3.0*C_DD) / 12.0) *
                (dpsifdS(l2, 0) / interpolated_grad_r[0]) *
                (dpsifdS(l, 0) / interpolated_grad_r[0]) * W * J;
            }

            // Loop over the velocity components
            for (unsigned i2 = 0; i2 < n_dim; i2++)
            {
              // Get the unknown u index
              local_unknown = this->nodal_local_eqn(l2, u_index[i2]);


              // If not a boundary condition
              if (local_unknown >= 0)
              {
                if (i2 == 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    -(GKE_Phi_C * interpolated_x[1]*(1.0+C_DD) / 2.0) * psif(l2) *
                    (dpsifdS(l, 0) / interpolated_grad_r[0]) * W * J;
                  jacobian(local_eqn, local_unknown) +=
                    -interpolated_dhdr * psif(l2) * psif(l) * W * J;
                }

                if (i2 == 1)
                {
                  jacobian(local_eqn, local_unknown) +=
                    psif(l2) * psif(l) * W * J;
                }
              }
            }
          }
        }
      }

      if (film_flag == true)
      {
        // Loop over velocity components
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Read out the apprporiate local equation for the u residual and
          // Jacobian terms
          local_eqn = this->nodal_local_eqn(l, u_index[i]);
          // If not a boundary condition
          if (local_eqn >= 0)
          {
            // Add Residual terms
            residuals[local_eqn] +=

              (-Vis_r * interpolated_Pv +
               Ha /
                 ((1.0+7.0*C_DD)*interpolated_x[1] * interpolated_x[1] * interpolated_x[1])) *
              interpolated_n[i] * psif(l) * W * J;
            if (i == 0)
            {
              residuals[local_eqn] +=
                -Vis_r *
                ((interpolated_x[1]*(1.0+C_DD) / (2.0 * GKE_Tau_P)) *
                   interpolated_grad_Pv[0] +
                 interpolated_u[0]*(1.0-C_DD) / (interpolated_x[1] * GKE_Tau_C)) *
                psif(l) * W * J;
            }
            // We also need to worry about the jacobian terms
            if (flag)
            {
              // Loop over the nodes again, this time for the unknowns
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Get the unknown c_index
                local_unknown = this->nodal_local_eqn(l2, this->Pv_index[l2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    Vis_r * (-psif(l2) * interpolated_n[i]) * psif(l) * W * J;
                  if (i == 0)
                  {
                    jacobian(local_eqn, local_unknown) +=
                      -Vis_r *
                      ((interpolated_x[1]*(1.0+C_DD) / (2.0 * GKE_Tau_P)) *
                       dpsifdS(l2, 0)) *
                      psif(l) * W * J;
                  }
                }

                // Loop over the velocity components
                for (unsigned i2 = 0; i2 < n_dim; i2++)
                {
                  // Get the unknown u index
                  local_unknown = this->nodal_local_eqn(l2, u_index[i2]);


                  // If not a boundary condition
                  if (local_unknown >= 0)
                  {
                    if (i == 0)
                    {
                      if (i2 == 0)
                      {
                        jacobian(local_eqn, local_unknown) +=
                          -Vis_r *
                          (psif(l2)*(1.0-C_DD) / (interpolated_x[1] * GKE_Tau_C)) *
                          psif(l) * W * J;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Read out the apprporiate local equation for the arc length residual and
      // Jacobian terms
      local_eqn = this->nodal_local_eqn(l, this->Al_index[l]);

      if (local_eqn >= 0)
      {
        for (unsigned i = 0; i < n_dim; i++)
        {
          residuals[local_eqn] +=
            interpolated_grad_Al[i] * dpsifdS(l, i) * W * J;
        }
        if (flag)
        {
          // Loop over the nodes again, this time for the unknowns
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            // Get the unknown Al_index
            local_unknown = this->nodal_local_eqn(l2, this->Al_index[l2]);
            if (local_unknown >= 0)
            {
              for (unsigned i = 0; i < n_dim; i++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dpsifdS(l2, i) * dpsifdS(l, i) * W * J;
              }
            }
          }
        }
      }
    }
  }

  //=======================================================
  /// Overload the output function
  //=======================================================
  void InterfaceLubricationElement::output(std::ostream& outfile,
                                           const unsigned& n_plot)
  {
    outfile.precision(16);

    const unsigned el_dim = this->dim();
    const unsigned n_dim = this->nodal_dimension();
    const unsigned n_velocity = this->U_index_interface.size();

    // Set output Vector
    Vector<double> s(el_dim);
    Vector<double> n(n_dim);
    Vector<double> u(n_velocity);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, n_plot, s);
      // Get the outer unit normal
      this->outer_unit_normal(s, n);

      // double u_n = 0.0;
      for (unsigned i = 0; i < n_velocity; i++)
      {
        u[i] = this->interpolated_u(s, i);
      }
      // Output the x,y,u,v
      for (unsigned i = 0; i < n_dim; i++)
        outfile << this->interpolated_x(s, i) << " ";
      for (unsigned i = 0; i < n_dim; i++) outfile << u[i] << " ";

      // Output the vapour pressure
      outfile << interpolated_Pv(s) << " ";
      outfile << interpolated_Al(s) << " ";

      outfile << std::endl;
    }
  }
} // namespace oomph
