#include "simple_output_elements.h"

namespace oomph
{

 void CleanOutputAxisymmetricTTaylorHoodElement::output(std::ostream& outfile, const unsigned& n_plot)
  {
  // Vector of local coordinates
    Vector<double> s(2);

     outfile.precision(16);


    // Loop over plot points
    unsigned num_plot_points = nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, n_plot, s);

      // Coordinates
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < 2; i++)
      {
        outfile << interpolated_u_axi_nst(s, i) << " ";
      }

      // Pressure
      outfile << interpolated_p_axi_nst(s) << " ";

      outfile << std::endl;
    }
    
  }
} // namespace oomph
