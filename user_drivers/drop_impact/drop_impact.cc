// Driver for axisymmetric drop impact onto a vapour layer.
// This is adapted from various demo driver codes in oomph-lib

// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// Constitutive law headers
#include "constitutive.h"

// Solid headers
#include "solid.h"

// The mesh headers
#include "meshes/triangle_mesh.h"

// My vapour element headers
#include "drop_impact_elements.h"

using namespace std;

using namespace oomph;

//==start_of_namespace====================================================
/// Namespace for physical parameters
//========================================================================
namespace Global_Physical_Variables {
/// Reynolds number
double Re = 32.0;

/// Product of Reynolds number and inverse of Froude number (Gravity Number)
double ReInvFr = 0.0;

/// Capillary number
double Ca = 0.0674;

/// Viscosity Ratio
double Vis_R = 0.0016;

///  Dim-less Hamaker number
double A_h = 0.0;

// Gas atmospheric pressure
// For incompressible flow, this can be set to 0.
double vapour_pressure_0 = 0.0;

/// Initial position for drop centre
double Drop_start_height = 1.1;

/// Initial Drop Radius
double Drop_radius = 1.0;

/// Initial velocity in z direction, so should be negative for an impact.
double Initial_velocity = -1.0;

/// Maximum time
double t_max = 3.5;

/// Use Gas Kinetic Effects?
bool Use_GKE = false;

/// Non-dim Mean free path of gas
double Mean_Free_Path = 8.625e-5;

/// Drop Drop Simulation?
bool Drop_Drop = false;

// Never need to change the rest of these

/// Strouhal number
double St = 1.0;

/// Womersley number (Reynolds x Strouhal, computed automatically)
double ReSt;

/// Direction of gravity
Vector<double> G(3);

/// Pseudo-solid Poisson ratio
double Nu = 0.2;

// GKE fuctions from Chubynsky et al 2020 for alpha=1
void GKE_Tau_P_fct(const double &Kn, double &result) { result = 1.0; }

void GKE_Tau_C_fct(const double &Kn, double &result) {
  double delta = sqrt(MathematicalConstants::Pi) / (2.0 * Kn);
  result = 1.0 + 2.0 * Kn *
                     (1.1466 - 0.095 * exp(-0.747 * delta) -
                      0.0516 * exp(-3.724 * delta));
}

void GKE_Phi_P_fct(const double &Kn, double &result) {
  if (Kn > 0.0) {
    result = 1.0 + 6.88 * Kn +
             (6.0 * Kn / MathematicalConstants::Pi) *
                 log1p(2.76 * Kn + 0.127 * Kn * Kn);
  } else {
    result = 1.0 - 6.88 * Kn -
             (6.0 * Kn / MathematicalConstants::Pi) *
                 log1p(-2.76 * Kn + 0.127 * Kn * Kn);
  }
}

void GKE_Phi_C_fct(const double &Kn, double &result) { result = 1.0; }

// Alternative GKE function from Zhang and Law 2011
void GKE_Phi_P_dd_fct(const double &Kn, double &result) {
  if (Kn < 1.0) {
    result = 1.0 + 6.0966 * Kn + 0.9650 * Kn * Kn + 0.6967 * Kn * Kn * Kn;
  } else {
    result = 8.7583 * pow(Kn, 1.1551);
  }
}

} // namespace Global_Physical_Variables

//==start_of_namespace====================================================
/// Namespace for Simulation Settings
//========================================================================
namespace Global_Sim_Settings {
/// Result Folder
string result_folder = "RESLT";

// Minimum allowed height, simulation stopped if reached.
double Min_height_cutoff = 1.0e-8;

// Maximum timestep
double max_dt = 0.0005;

// Maximum timestep during impact
double max_dt_impact = 0.0001;

// Minimum timestep
double min_dt = 1.0e-12;

// target error for adaptive timestepping
double epsilon_t = 1.0e-5;

// Max element area for initial mesh
double initial_max_element_area = 0.01;

// Number of points on free surface in initial mesh
unsigned initial_n_surface_points = 100;

// Number of timesteps per adaption
unsigned steps_per_adapt = 100;

// Number of steps per full output - set very large to supress full outputs by default
unsigned steps_per_full_output = 100000000;

// Start of Impact (where max_dt_impact is used)
double impact_start_time = 0.09;

// End of impact
double impact_end_time = 0.17;

// Targets for spatial adaptivity
// Max error (elements refined if error is above this)
double max_permitted_error = 5.0e-5;
// Min error (elements unrefined if error is below this)
double min_permitted_error = 1.0e-5;
// Max element size
double max_element_size = 0.01;
// Min element size
double min_element_size = 1.0e-7;

// Tolerance for refining boundary based on curvature
double polyline_refinement_tolerance = 0.002;   // Default 0.08
// Tolerance for unrefining boundary based on curvature
double polyline_unrefinement_tolerance = 0.001; // Default 0.04

// Lubrication equation is solved  where 
// the z component of the normal to the drop is less than -film_cutoff
double film_cutoff = 0.1; 

} // End of namespace Global_Sim_Settings

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//==start_of_problem_class================================================
/// Impact of a drop onto a vapour film
//========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class ImpactProblem : public Problem {
public:
  /// Constructor
  ImpactProblem();

  /// Destructor
  ~ImpactProblem() {
    delete outer_boundary_circle_pt;

    for (unsigned i = 0; i < 2; i++) {
      delete boundary_polyline_pt[i];
    }
    boundary_polyline_pt.clear();
    delete boundary_polygon_pt;
    delete Bulk_mesh_pt;
    delete Surface_mesh_pt;
    delete Constitutive_law_pt;
    delete External_pressure_data_pt;
    delete error_estimator_pt;
    delete this->time_stepper_pt();
  }

  /// Set initial conditions
  void set_initial_condition();

  /// Document the solution
  void doc_solution(DocInfo &doc_info, bool full_output, bool post_adapt);

  // Function to find minimum height
  double min_height();

  /// Global error norm for adaptive time-stepping
  double global_temporal_error_norm();

  void update_pressure_pins() { // Repin vapour pressure boundary conditions

    // Determine number of 1D interface elements in mesh
    const unsigned n_interface_element = Surface_mesh_pt->nelement();

    for (unsigned e = 0; e < n_interface_element; e++) {
      // Upcast from GeneralisedElement to the present element
      INTERFACE_ELEMENT *el_pt =
          dynamic_cast<INTERFACE_ELEMENT *>(Surface_mesh_pt->element_pt(e));

      // Get the outer unit normal - this is hacky, copying the output function
      unsigned n_dim = el_pt->nodal_dimension();
      unsigned el_dim = el_pt->dim();
      Vector<double> n(n_dim);
      Vector<double> s(el_dim);
      // Get local coordinates of plot point
      el_pt->get_s_plot(1, 3, s);
      // Get the outer unit normal
      el_pt->outer_unit_normal(s, n);
      
      // Is this element part of the film?
      bool film_flag = false;
      if (n[1] < -Global_Sim_Settings::film_cutoff) {
        film_flag = true;
      }
      // Pin the gas pressure if not in the film
      const unsigned nnodes = el_pt->nnode();
      for (unsigned i = 0; i < nnodes; i++) {
        el_pt->unpin_vapour_pressure(i);
        if (film_flag == false) {
          el_pt->fix_vapour_pressure(
              i, Global_Physical_Variables::vapour_pressure_0);
        }
      }
    }
    // Reassign equation numbers due to the re-pinning
    oomph_info << "\nNumber of equations: " << assign_eqn_numbers()
               << std::endl;
  }

private:
  /// Actions required before solve step (empty)
  void actions_before_newton_solve() {}

  /// No actions required after solve step (empty)
  void actions_after_newton_solve() {}

  /// Actions before adapt: Wipe the mesh of boundary elements
  void actions_before_adapt() {
    // Kill the  elements and wipe surface mesh
    delete_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

  } // end of actions_before_adapt

  /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
  void actions_after_adapt() {
    // Reassign the lagrangian coordinates (an updated lagrangian approach)
    Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

    // Create the boundary elements
    create_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

    // Setup the problem again -- remember that fluid mesh has been
    // completely rebuilt and its element's don't have any
    // pointers to Re etc. yet
    complete_problem_setup();

  } // end of actions_after_adapt

  /// \short Actions before the timestep: For maximum stability, reset
  /// the current nodal positions to be the "stress-free" ones.
  void actions_before_implicit_timestep() {
    Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

  /// Create free surface elements
  void create_boundary_elements();

  /// Delete boundary elements and wipe the associated mesh
  void delete_boundary_elements() {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++) {
      // Kill surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();

  } // end of delete_boundary_elements

  void complete_problem_setup() {
    // ----------------------------------------------------------------
    // Complete the problem setup to make the elements fully functional
    // ----------------------------------------------------------------

    // Determine number of bulk elements in mesh
    const unsigned n_element_bulk = Bulk_mesh_pt->nelement();

    // Loop over the bulk elements
    for (unsigned e = 0; e < n_element_bulk; e++) {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *el_pt = dynamic_cast<ELEMENT *>(Bulk_mesh_pt->element_pt(e));

      // Set the Reynolds number
      el_pt->re_pt() = &Global_Physical_Variables::Re;

      // Set the Womersley number
      el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

      // Set the product of the Reynolds number and the inverse of the
      // Froude number
      el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

      // Set the direction of gravity
      el_pt->g_pt() = &Global_Physical_Variables::G;

      // Set the constitutive law
      el_pt->constitutive_law_pt() = Constitutive_law_pt;

    } // End of loop over bulk elements

    // Determine number of 1D interface elements in mesh
    const unsigned n_interface_element = Surface_mesh_pt->nelement();
    // Counters to locate the node at the top of the drop
    unsigned top_element = 0;
    unsigned top_node = 0;
    double max_z_found = -1000.0;
    unsigned bottom_element = 0;
    unsigned bottom_node = 0;
    double min_z_found = 1000.0;
    // Loop over the interface elements
    for (unsigned e = 0; e < n_interface_element; e++) {
      // Upcast from GeneralisedElement to the present element
      INTERFACE_ELEMENT *el_pt =
          dynamic_cast<INTERFACE_ELEMENT *>(Surface_mesh_pt->element_pt(e));

      // Set the Strouhal number
      el_pt->st_pt() = &Global_Physical_Variables::St;

      // Set the Capillary number
      el_pt->ca_pt() = &Global_Physical_Variables::Ca;

      // Pass the Data item that contains the single external pressure value
      el_pt->set_external_pressure_data(External_pressure_data_pt);

      // Set the viscosity ratio
      el_pt->vapour_viscosity_ratio_pt() = &Global_Physical_Variables::Vis_R;

      // Set the Hamaker number
      el_pt->hamaker_pt() = &Global_Physical_Variables::A_h;

      // Set the MFP
      el_pt->mean_free_path_pt() = &Global_Physical_Variables::Mean_Free_Path;

      // Set the film cutoff
      el_pt->film_cutoff_pt() = &Global_Sim_Settings::film_cutoff;

      // Set the drop-drop flag
      el_pt->drop_drop_flag_pt() = &Global_Physical_Variables::Drop_Drop;

      // Set the GKE factors
      if (Global_Physical_Variables::Use_GKE) {
        if (Global_Physical_Variables::Drop_Drop == false) {
          el_pt->gke_Phi_P_fct_pt() = &Global_Physical_Variables::GKE_Phi_P_fct;
          el_pt->gke_Phi_C_fct_pt() = &Global_Physical_Variables::GKE_Phi_C_fct;
          el_pt->gke_Tau_P_fct_pt() = &Global_Physical_Variables::GKE_Tau_P_fct;
          el_pt->gke_Tau_C_fct_pt() = &Global_Physical_Variables::GKE_Tau_C_fct;
        }
        if (Global_Physical_Variables::Drop_Drop == true) {
          el_pt->gke_Phi_P_fct_pt() = &Global_Physical_Variables::GKE_Phi_P_fct;
          el_pt->gke_Phi_C_fct_pt() = &Global_Physical_Variables::GKE_Phi_C_fct;
          el_pt->gke_Tau_P_fct_pt() = &Global_Physical_Variables::GKE_Tau_P_fct;
          el_pt->gke_Tau_C_fct_pt() = &Global_Physical_Variables::GKE_Tau_P_fct;
        }
      }

      // Loop over nodes
      const unsigned nnodes = el_pt->nnode();
      for (unsigned i = 0; i < nnodes; i++) {
        Node *nod_pt = el_pt->node_pt(i);
        double r = nod_pt->x(0);
        double z = nod_pt->x(1);
        // Identify if node is at the top of the drop
        if (r < 1.0e-15) {
          if (z > max_z_found) {
            // Save number of identified node
            max_z_found = z;
            top_element = e;
            top_node = i;
          }

          if (z < min_z_found) {
            // Save number of identified node
            min_z_found = z;
            bottom_element = e;
            bottom_node = i;
          }
        }
      }

    } // End of loop over interface elements

    // Pin identified Vapour pressure in identified node
    // Upcast from GeneralisedElement to the present element
    INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT *>(
        Surface_mesh_pt->element_pt(top_element));
    Node *nod_pt = el_pt->node_pt(top_node);
    double r = nod_pt->x(0);
    double z = nod_pt->x(1);
    // Fix vapour pressure
    el_pt->fix_vapour_pressure(top_node,
                               Global_Physical_Variables::vapour_pressure_0);
    // Fix the arc length
    el_pt->fix_arc_length(top_node, 1.0);
    // Document to check this has happened
    oomph_info << "Pinned vapour="
               << Global_Physical_Variables::vapour_pressure_0
               << "and AL=1 at r=" << r << " and z=" << z << std::endl;

    // Pin Al at bottom in identified node
    // Upcast from GeneralisedElement to the present element
    el_pt = dynamic_cast<INTERFACE_ELEMENT *>(
        Surface_mesh_pt->element_pt(bottom_element));
    nod_pt = el_pt->node_pt(bottom_node);
    r = nod_pt->x(0);
    z = nod_pt->x(1);
    // Fix the arc length
    el_pt->fix_arc_length(bottom_node, 0.0);
    // Document to check this has happened
    oomph_info << "Pinned AL=0 at r=" << r << " and z=" << z << std::endl;

    // --------------------------------------------
    // Set the boundary conditions for this problem
    // --------------------------------------------

    // All nodes are free by default -- just pin the ones that have
    // Dirichlet conditions here

    // Fix Symmetry boundary
    // Find number of nodes on boundary
    const unsigned n_boundary_node =
        Bulk_mesh_pt->nboundary_node(fixed_boundary_id);
    // Loop over nodes on fixed boundary
    for (unsigned n = 0; n < n_boundary_node; n++) {
      // Pin Value of radial velocity to 0
      Bulk_mesh_pt->boundary_node_pt(fixed_boundary_id, n)->pin(0);

      // Set value to 0
      Bulk_mesh_pt->boundary_node_pt(fixed_boundary_id, n)->set_value(0, 0.0);
      // Pin radial position
      Bulk_mesh_pt->boundary_node_pt(fixed_boundary_id, n)->pin_position(0);
    }

    // Pin all azimuthal velocities
    const unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++) {
      Bulk_mesh_pt->node_pt(n)->pin(2);
    }
  }

  /// Pointer to the (specific) "bulk" mesh
  RefineableSolidTriangleMesh<ELEMENT> *Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh *Surface_mesh_pt;

  // Pointer to the external pressure
  Data *External_pressure_data_pt;

  // Pointer to the constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt;

  // Pointers for mesh construction
  Circle *outer_boundary_circle_pt;
  Vector<TriangleMeshCurveSection *> boundary_polyline_pt;
  TriangleMeshPolygon *boundary_polygon_pt;
  TriangleMeshClosedCurve *closed_curve_pt;

  // Pointer to error estimator
  Z2ErrorEstimator *error_estimator_pt;

  // ids of the parts of the boundary
  // Free boundary at surface of drop
  const unsigned free_boundary_id = 0;
  // Flat symmetry boundary at r=0
  const unsigned fixed_boundary_id = 1;

}; // End of problem class

//==start_of_constructor==================================================
/// Constructor for single fluid free surface problem
//========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::ImpactProblem() {
  // Allocate the timestepper (this constructs the time object as well)
  add_time_stepper_pt(new BDF<2>(true));
#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  linear_solver_pt() = new MumpsSolver;
#endif
  // Don't reject timesteps above tolerance
  Keep_temporal_error_below_tolerance = false;

  // Coordinate along surface
  Vector<double> zeta(1);

  // Position vector on GeomObject
  Vector<double> posn(2);

  // Circle object used to define the outer boundary
  outer_boundary_circle_pt =
      new Circle(0.0, Global_Physical_Variables::Drop_start_height,
                 Global_Physical_Variables::Drop_radius);
  // Pointer to the closed curve that defines the outer boundary
  closed_curve_pt = 0;

  // Provide storage for pointers to the two parts of the curvilinear boundary
  boundary_polyline_pt.resize(2);

  // Free boundary
  // Define start and end angle, anticlockwise
  double zeta_start = -0.5 * MathematicalConstants::Pi;
  double zeta_end = 0.5 * MathematicalConstants::Pi;
  // Contruct the points on the free boundary
  unsigned npoints = Global_Sim_Settings::initial_n_surface_points;
  double unit_zeta = MathematicalConstants::Pi / double(npoints - 1);
  Vector<Vector<double>> free_boundary_vertex(npoints);
  for (unsigned ipoint = 0; ipoint < npoints; ipoint++) {
    // Resize the vector
    free_boundary_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = zeta_start + unit_zeta * double(ipoint);
    outer_boundary_circle_pt->position(zeta, posn);
    free_boundary_vertex[ipoint] = posn;
  }
  boundary_polyline_pt[0] =
      new TriangleMeshPolyLine(free_boundary_vertex, free_boundary_id);

  // Flat boundary
  unsigned n_poly_points = 2;
  Vector<Vector<double>> vertex_coord(n_poly_points);
  for (unsigned i = 0; i < n_poly_points; i++) {
    vertex_coord[i].resize(2);
  }
  zeta[0] = zeta_end;
  outer_boundary_circle_pt->position(zeta, posn);
  vertex_coord[0][0] = posn[0];
  vertex_coord[0][1] = posn[1];
  zeta[0] = zeta_start;
  outer_boundary_circle_pt->position(zeta, posn);
  vertex_coord[1][0] = posn[0];
  vertex_coord[1][1] = posn[1];

  // Build the 1st boundary polyline
  boundary_polyline_pt[1] =
      new TriangleMeshPolyLine(vertex_coord, fixed_boundary_id);

  // Combine to curvilinear boundary
  //--------------------------------
  boundary_polygon_pt = new TriangleMeshPolygon(boundary_polyline_pt);
  // Set the refinement tolerances
  boundary_polygon_pt->set_polyline_refinement_tolerance(
      Global_Sim_Settings::polyline_refinement_tolerance);

  boundary_polygon_pt->set_polyline_unrefinement_tolerance(
      Global_Sim_Settings::polyline_unrefinement_tolerance);

  // Set the pointer
  closed_curve_pt = boundary_polygon_pt;

  // Now build the mesh
  //===================

  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters
  TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

  // Specify the maximum area element
  triangle_mesh_parameters.element_area() =
      Global_Sim_Settings::initial_max_element_area;

  // Store as the problem's main "Bulk" mesh
  Bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Bulk_mesh_pt->max_permitted_error() =
      Global_Sim_Settings::max_permitted_error;
  Bulk_mesh_pt->min_permitted_error() =
      Global_Sim_Settings::min_permitted_error;
  Bulk_mesh_pt->max_element_size() = Global_Sim_Settings::max_element_size;
  Bulk_mesh_pt->min_element_size() = Global_Sim_Settings::min_element_size;
  Bulk_mesh_pt->set_print_level_timings_adaptation(3);

  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Surface_mesh_pt = new Mesh;

  // -----------------------------
  // Create the interface elements
  // -----------------------------

  create_boundary_elements();

  // Add the two sub-meshes to the problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  // Define a constitutive law for the solid equations: generalised Hookean
  Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);

  // Create a Data object whose single value stores the dummy external pressure
  External_pressure_data_pt = new Data(1);
  External_pressure_data_pt->pin(0);
  External_pressure_data_pt->set_value(0, 0.0);

  complete_problem_setup();

  // Setup equation numbering scheme
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor

//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
double ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::global_temporal_error_norm() {
  double global_position_error = 0.0;

  double global_velocity_error = 0.0;

  // Find out how many nodes there are in the problem
  unsigned n_node = mesh_pt()->nnode();

  // Loop over the nodes and calculate the estimated error in the values
  for (unsigned i = 0; i < n_node; i++) {
    // Get error in solution: Difference between predicted and actual
    // value for position value 0
    double error = mesh_pt()
                       ->node_pt(i)
                       ->position_time_stepper_pt()
                       ->temporal_error_in_position(mesh_pt()->node_pt(i), 0);

    // Add the square of the individual error to the global error
    global_position_error += error * error;

    // Get error in solution: Difference between predicted and actual
    // value for position value 1
    error = mesh_pt()
                ->node_pt(i)
                ->position_time_stepper_pt()
                ->temporal_error_in_position(mesh_pt()->node_pt(i), 1);

    // Add the square of the individual error to the global error
    global_position_error += error * error;

    error = mesh_pt()->node_pt(i)->time_stepper_pt()->temporal_error_in_value(
        mesh_pt()->node_pt(i), 0);

    global_velocity_error += error * error;

    error = mesh_pt()->node_pt(i)->time_stepper_pt()->temporal_error_in_value(
        mesh_pt()->node_pt(i), 1);

    global_velocity_error += error * error;
  }

  // Divide by the number of nodes
  global_position_error /= 2.0 * double(n_node);

  global_velocity_error /= 2.0 * double(n_node);

  // Return sum of square roots of errors
  return sqrt(global_position_error) + sqrt(global_velocity_error);
} // end of global_temporal_error_norm

//============start_of_create_boundary_elements===============
/// Create boundary elements
//=======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::create_boundary_elements() {
  unsigned n_element = Bulk_mesh_pt->nboundary_element(0);
  // Loop over those elements adjacent to the free surface
  for (unsigned e = 0; e < n_element; e++) {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT *bulk_element_pt =
        dynamic_cast<ELEMENT *>(Bulk_mesh_pt->boundary_element_pt(0, e));

    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(0, e);

    // Create the interface element
    INTERFACE_ELEMENT *interface_element_pt =
        new INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(0);
  }

} // end of create_boundary_elements

//==start_of_set_initial_condition========================================
/// \short Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::set_initial_condition() {
  // Determine number of nodes in mesh
  const unsigned n_node = mesh_pt()->nnode();

  // Loop over all nodes in mesh
  for (unsigned n = 0; n < n_node; n++) {
    // Loop over the velocity components
    // Assign radial and angular to 0, set vertical to Initial velocity
    mesh_pt()->node_pt(n)->set_value(0, 0.0);
    mesh_pt()->node_pt(n)->set_value(
        1, Global_Physical_Variables::Initial_velocity);
    mesh_pt()->node_pt(n)->set_value(2, 0.0);
  }

  // Loop over surface mesh

  const unsigned n_interface_element = Surface_mesh_pt->nelement();
  for (unsigned e = 0; e < n_interface_element; e++) {
    INTERFACE_ELEMENT *el_pt =
        dynamic_cast<INTERFACE_ELEMENT *>(Surface_mesh_pt->element_pt(e));

    const unsigned nnodes = el_pt->nnode();
    for (unsigned i = 0; i < nnodes; i++) {
      el_pt->set_vapour_pressure(i,
                                 Global_Physical_Variables::vapour_pressure_0);
    }
  }

} // End of set_initial_condition

// Compute min height
template <class ELEMENT, class INTERFACE_ELEMENT>
double ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::min_height() {
  double min_value = 1000.0;
  unsigned n_plot = 20;
  // Loop over surface elemnets
  const unsigned n_interface_element = Surface_mesh_pt->nelement();
  for (unsigned e = 0; e < n_interface_element; e++) {
    INTERFACE_ELEMENT *el_pt =
        dynamic_cast<INTERFACE_ELEMENT *>(Surface_mesh_pt->element_pt(e));
    const unsigned el_dim = el_pt->dim();
    Vector<double> s(el_dim);
    // Loop over points on surface element
    unsigned num_plot_points = el_pt->nplot_points(n_plot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++) {
      // Get local coordinates of plot point
      el_pt->get_s_plot(iplot, n_plot, s);
      // Check if z is less than lowest previously found
      double z = el_pt->interpolated_x(s, 1);
      if (z < min_value) {
        min_value = z;
      }
    }
  }
  return min_value;
} // End of min_height

//==start_of_doc_solution=================================================
/// Document the solution
//========================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, INTERFACE_ELEMENT>::doc_solution(DocInfo &doc_info,
                                                             bool full_output,
                                                             bool post_adapt) {
  // Output the time
  oomph_info << "Outputting at time " << time_pt()->time() << std::endl;

  ofstream some_file;
  char filename_surface[100];
  char filename_bulk[100];
  char filename_size[100];
  if (post_adapt) {
    sprintf(filename_bulk, "%s/elements_%i_post_adapt.dat",
            doc_info.directory().c_str(), doc_info.number());

    sprintf(filename_size, "%s/elements_size_%i_post_adapt.dat",
            doc_info.directory().c_str(), doc_info.number());

    sprintf(filename_surface, "%s/surface_elements_%i_post_adapt.dat",
            doc_info.directory().c_str(), doc_info.number());
  } else {
    sprintf(filename_bulk, "%s/elements_%i.dat", doc_info.directory().c_str(),
            doc_info.number());

    sprintf(filename_size, "%s/elements_size_%i.dat",
            doc_info.directory().c_str(), doc_info.number());

    sprintf(filename_surface, "%s/surface_elements_%i.dat",
            doc_info.directory().c_str(), doc_info.number());
  }
  if (full_output) {
    // Output bulk elements
    some_file.open(filename_bulk);
    Bulk_mesh_pt->output(some_file, 3);
    some_file.close();

    // Output element areas
    some_file.open(filename_size);
    int n_elements = this->Bulk_mesh_pt->nelement();
    for (int j = 0; j < n_elements; j++) {
      ELEMENT *el_pt =
          dynamic_cast<ELEMENT *>(this->Bulk_mesh_pt->element_pt(j));
      some_file << j << " " << el_pt->size() << std::endl;
    }
    some_file.close();
  }

  // Save the vapour boundary
  some_file.open(filename_surface);
  Surface_mesh_pt->output(some_file, 3);
  some_file.close();
} // End of doc_solution

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//==start_of_main=========================================================
/// Driver code for axisym drop impact with vapour layer
//========================================================================
int main(int argc, char *argv[]) {
#ifdef OOMPH_HAS_MPI

  // Initialise MPI
  MPI_Helpers::init(argc, argv);

  // oomph_mpi_output.restrict_output_to_single_processor();
#endif

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
      "--folder", &Global_Sim_Settings::result_folder, "Result Folder");

  CommandLineArgs::specify_command_line_flag(
      "--reynolds", &Global_Physical_Variables::Re, "Reynolds Number");

  CommandLineArgs::specify_command_line_flag(
      "--gravity", &Global_Physical_Variables::ReInvFr, "Gravity Number");

  CommandLineArgs::specify_command_line_flag(
      "--capilliary", &Global_Physical_Variables::Ca, "Capilliary Number");

  CommandLineArgs::specify_command_line_flag(
      "--viscosity", &Global_Physical_Variables::Vis_R, "Viscosity Ratio");

  CommandLineArgs::specify_command_line_flag(
      "--hamaker", &Global_Physical_Variables::A_h, "Dim-less Hamaker number");

  CommandLineArgs::specify_command_line_flag(
      "--mfp", &Global_Physical_Variables::Mean_Free_Path,
      "Dim-less Mean-free-path");

  CommandLineArgs::specify_command_line_flag("--usegke", "Use GKE functions");

  CommandLineArgs::specify_command_line_flag("--dropdrop", "Drop Drop Impact");

  CommandLineArgs::specify_command_line_flag(
      "--tmax", &Global_Physical_Variables::t_max, "Max time");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // Finalise parameters
  // Compute the Womersley number
  Global_Physical_Variables::ReSt =
      Global_Physical_Variables::Re * Global_Physical_Variables::St;
  // Set direction of gravity (vertically downwards)
  Global_Physical_Variables::G[0] = 0.0;
  Global_Physical_Variables::G[1] = -1.0;
  Global_Physical_Variables::G[2] = 0.0;

  if (CommandLineArgs::command_line_flag_has_been_set("--usegke")) {
    Global_Physical_Variables::Use_GKE = true;
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--dropdrop")) {
    Global_Physical_Variables::Drop_Drop = true;
  }

  oomph_info << "Driver: Constructing Problem" << std::endl;
  // Create problem with Axisym Taylor hood element, with pseudosolid mesh
  // update and projectability for remeshing
  ImpactProblem<
      ProjectableAxisymmetricTaylorHoodElement<PseudoSolidNodeUpdateElement<
          CleanOutputAxisymmetricTTaylorHoodElement, TPVDElement<2, 3>>>,
      ElasticAxisymmetricInterfaceLubricationElement<
          ProjectableAxisymmetricTaylorHoodElement<PseudoSolidNodeUpdateElement<
              CleanOutputAxisymmetricTTaylorHoodElement, TPVDElement<2, 3>>>>>
      problem;

  // Initialise DocInfo object
  DocInfo doc_info;

  // Set output directory
  doc_info.set_directory(Global_Sim_Settings::result_folder);

  // Initialise counter for solutions
  doc_info.number() = 0;

  // Open a trace file
  ofstream trace_file;
  char filename[100];
  sprintf(filename, "%s/trace.dat", doc_info.directory().c_str());
  trace_file.open(filename);
  trace_file.precision(16);

  // Print Trace file header
  trace_file << "---------------------------------------------" << std::endl;
  trace_file << "----------------- Drop Impact ---------------" << std::endl;
  trace_file << "---------------------------------------------" << std::endl;
  trace_file << "Output Folder:       " << Global_Sim_Settings::result_folder
             << std::endl;
  trace_file << "Reynolds Number:     " << Global_Physical_Variables::Re
             << std::endl;
  trace_file << "Gravity Number:      " << Global_Physical_Variables::ReInvFr
             << std::endl;
  trace_file << "Capilliary Number:   " << Global_Physical_Variables::Ca
             << std::endl;
  trace_file << "Viscosity Ratio:     " << Global_Physical_Variables::Vis_R
             << std::endl;
  trace_file << "NDim Hamaker Number: " << Global_Physical_Variables::A_h
             << std::endl;
  trace_file << "Centre Start Height: "
             << Global_Physical_Variables::Drop_start_height << std::endl;
  trace_file << "Initial Drop Radius: "
             << Global_Physical_Variables::Drop_radius << std::endl;
  trace_file << "Initial Velocity:    "
             << Global_Physical_Variables::Initial_velocity << std::endl;
  trace_file << "Maximum time:        " << Global_Physical_Variables::t_max
             << std::endl;
  trace_file << "Using GKE:           " << Global_Physical_Variables::Use_GKE
             << std::endl;
  trace_file << "Ndim Mean Free Path: "
             << Global_Physical_Variables::Mean_Free_Path << std::endl;
  trace_file << "Drop-Drop:           " << Global_Physical_Variables::Drop_Drop
             << std::endl;
  trace_file << "Min allowed height:  "
             << Global_Sim_Settings::Min_height_cutoff << std::endl;
  trace_file << "Max Timestep:        " << Global_Sim_Settings::max_dt
             << std::endl;
  trace_file << "Max Timestep Impact: " << Global_Sim_Settings::max_dt_impact
             << std::endl;
  trace_file << "Min Timestep:        " << Global_Sim_Settings::min_dt
             << std::endl;
  trace_file << "Target Time Error:   " << Global_Sim_Settings::epsilon_t
             << std::endl;
  trace_file << "Start Element Area:  "
             << Global_Sim_Settings::initial_max_element_area << std::endl;
  trace_file << "Init Surface points: "
             << Global_Sim_Settings::initial_n_surface_points << std::endl;
  trace_file << "Steps per adapt:     " << Global_Sim_Settings::steps_per_adapt
             << std::endl;
  trace_file << "Steps per full out:  "
             << Global_Sim_Settings::steps_per_full_output << std::endl;
  trace_file << "Max mesh error:      "
             << Global_Sim_Settings::max_permitted_error << std::endl;
  trace_file << "Min mesh error:      "
             << Global_Sim_Settings::min_permitted_error << std::endl;
  trace_file << "max element size:    " << Global_Sim_Settings::max_element_size
             << std::endl;
  trace_file << "min element size:    " << Global_Sim_Settings::min_element_size
             << std::endl;
  trace_file << "bound refine tol:    "
             << Global_Sim_Settings::polyline_refinement_tolerance << std::endl;
  trace_file << "bound unrefine tol:  "
             << Global_Sim_Settings::polyline_unrefinement_tolerance
             << std::endl;
  trace_file << "Impact Start Time:   "
             << Global_Sim_Settings::impact_start_time << std::endl;
  trace_file << "Impact End Time:     " << Global_Sim_Settings::impact_end_time
             << std::endl;
  trace_file << "---------------------------------------------" << std::endl;

  // Headers for trace file
  trace_file << "VARIABLES=\"timestep\",\"time\",\"dt\",\"Global Error\",\"min "
                "z\",\"N Elements\""
             << std::endl;

  // Set initial conditions
  problem.set_initial_condition();

  // variables for timestepper
  double min_height;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(Global_Sim_Settings::max_dt);

  problem.maximum_dt() = Global_Sim_Settings::max_dt;
  problem.minimum_dt() = Global_Sim_Settings::min_dt;

  // Initial timestep: Use the one used when setting up the initial
  // condition
  double dt = problem.time_pt()->dt();
  // Initialise the previous velocity values and nodal positions
  // for timestepping corresponding to an impulsive start
  problem.assign_initial_values_impulsive();

  // Could output initial condition here. Instead just increment the info number
  doc_info.number()++;
  // Run the unsteady simulation

  // Counters for mesh adaption and output
  unsigned adaptcount = 1;
  bool before_impact = true;
  bool before_end_impact = true;
  bool adapt_flag = false;
  unsigned full_out_count = 1;
  // Loop until final time
  while (problem.time_pt()->time() < Global_Physical_Variables::t_max) {
    oomph_info << "Driver: Computing Step " << doc_info.number() << " at time "
               << problem.time_pt()->time() << std::endl;

    // Adapt periodically
    if (adaptcount == Global_Sim_Settings::steps_per_adapt) {
      adapt_flag = true;
    }
    if (before_impact) {
      if (problem.time_pt()->time() > Global_Sim_Settings::impact_start_time) {
        before_impact = false;
        problem.maximum_dt() = Global_Sim_Settings::max_dt_impact;
      }
    } else if (before_end_impact) {
      if (problem.time_pt()->time() > Global_Sim_Settings::impact_end_time) {
        before_end_impact = false;
        problem.maximum_dt() = Global_Sim_Settings::max_dt;
      }
    }

    if (adapt_flag == true) {
      adapt_flag = false;
      adaptcount = 0;
      oomph_info << "Driver: Adapting:" << std::endl;
      problem.adapt();
      oomph_info << "Driver: Finished Adaption" << std::endl;
      // Uncomment to output full solution after adaption (for debugging)
      //      problem.doc_solution(doc_info, true, true);
    }
    adaptcount = adaptcount + 1;

    problem.update_pressure_pins();
    // Take one adaptive timestep
    double dt_next = problem.adaptive_unsteady_newton_solve(
        dt, Global_Sim_Settings::epsilon_t);

    oomph_info << "Driver: Finished Step" << std::endl;

    // Doc solution
    oomph_info << "Driver: Docing Solution" << std::endl;
    if (full_out_count == Global_Sim_Settings::steps_per_full_output) {
      problem.doc_solution(doc_info, true, false);
      full_out_count = 0;
    } else {
      problem.doc_solution(doc_info, false, false);
    }
    full_out_count = full_out_count + 1;

    oomph_info << "Driver: Trace" << std::endl;

    trace_file << doc_info.number() << " " << problem.time_pt()->time() << " "
               << problem.time_pt()->dt() << " "
               << problem.global_temporal_error_norm() << " "
               << problem.min_height() << " " << problem.mesh_pt()->nelement()
               << std::endl;

    // Increment counter for solutions
    doc_info.number()++;
    dt = dt_next;

    min_height = problem.min_height();
    // Exit loop if contact found
    if (min_height < Global_Sim_Settings::Min_height_cutoff) {
      oomph_info << "Driver: Contact detected, ending simulation. Min height="
                 << min_height << std::endl;
      break;
    }
    // Exit if bounce detected
    if (min_height > 1.1 * (Global_Physical_Variables::Drop_start_height -
                            Global_Physical_Variables::Drop_radius)) {
      oomph_info << "Driver: Rebound detected, ending simulation. Min height="
                 << min_height << std::endl;
      break;
    }
    oomph_info << "Driver: Finished Timestep Loop" << std::endl;
  } // End of timestepping loop

  // Close trace file
  trace_file.close();

  // Finalise MPI
#ifdef OOMPH_HAS_MPI

  MPI_Helpers::finalize();

#endif
} // End of main
