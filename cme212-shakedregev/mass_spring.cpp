/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"
#include "SpaceSearcher.hpp"



// Gravity in meters/sec^2
static constexpr double grav = 9.81;
double damping;
double L;
double K;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double K, L;
  EdgeData() : K(100), L(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

struct reset_position
{
  reset_position(const double dt) : dt_(dt) {}
  void operator()(Node n) {
    /* Prevent the cloth from simply falling to infinity, we can constrain two corners */
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1 , 0 , 0)) {
      n.value().vel = Point(0);
    }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    else{ //not needed, but to stress nothing changes
      n.position() += n.value().vel * dt_;
    }
  }
  const double dt_;
};

template <typename F>
struct reset_velocity
{
  reset_velocity(const double t, const double dt, F force) : t_(t), dt_(dt), force_(force) {}
  // Compute the t+dt velocity
  void operator()(Node n) {
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m because F=ma
    n.value().vel += force_(n, t_) * dt_ / n.value().mass;
  }
  const double t_;
  const double dt_;
  F force_;
};


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
// template <typename G, typename F>
// double symp_euler_step(G& g, double t, double dt, F force) {
//   // Compute the t+dt position
//   for (auto it = g.node_begin(); it != g.node_end(); ++it) {
//     auto n = *it;

//     // Update the position of the node according to its velocity
//     // x^{n+1} = x^{n} + v^{n} * dt
//     n.position() += n.value().vel * dt;
//   }

//   // Compute the t+dt velocity
//   for (auto it = g.node_begin(); it != g.node_end(); ++it) {
//     auto n = *it;

//     // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
//     n.value().vel += force(n, t) * (dt / n.value().mass);
//   }

//   return t + dt;
// }
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
//parallel code
  thrust::for_each(thrust::omp::par,g.node_begin(),g.node_end(),reset_position(dt));
  // /* Apply constraints */
  constraints(g, t);
  thrust::for_each(thrust::omp::par,g.node_begin(),g.node_end(),reset_velocity<F>(t,dt,force));
  return t + dt;

//serial code
  // for (auto it = g.node_begin(); it != g.node_end(); ++it) {
  //   auto n = *it;

  //   /* Prevent the cloth from simply falling to infinity, we can constrain two corner. */
  //   if (n.position() == Point(0, 0, 0) || n.position() == Point(1 , 0 , 0)) {
  //     n.value().vel = Point(0);
  //   }

  //   // Update the position of the node according to its velocity
  //   // x^{n+1} = x^{n} + v^{n} * dt
  //   else //this is not strictly needed because we would just multiply by 0, clarifying position doesn't change
  //   {
  //     n.position() += n.value().vel * dt;
  //   }
  // }
  // /* Apply constraints */
  // constraints(g, t);

  // // Compute the t+dt velocity
  // for (auto it = g.node_begin(); it != g.node_end(); ++it) {
  //   auto n = *it;

  //   // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
  //   n.value().vel += force(n, t) * (dt / n.value().mass);
  // }
  // return t + dt;

}
/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    // return Point(0);
    /* Preventing the cloth from simply falling to infinity, we constrain two corners. */
    Point force = Point(0, 0, 0);
    if (n.position () ==  Point (0,0,0) || n.position () ==  Point (1,0,0))
      return  force;
    force-= n.value().mass*Point(0, 0,grav);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point x_j = (*it).node2().position();
      force -= (*it).value().K * ((n.position() - x_j) / norm(n.position() - x_j)) * (norm(n.position() - x_j) - (*it).value().L);
      // force -= K * ((n.position() - x_j) / norm(n.position() - x_j)) * (norm(n.position() - x_j) - L);
    }
    return force;
  }
};

struct GravityForce {
  /** Return the force of gravity applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; // Quiet compiler warning
    return Point(0, 0, -n.value().mass*grav);
  }
};

/** Force function object that implements the spring force for HW2 #4. */
struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; // Quiet compiler warning
    Point force_spring = Point(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) 
    {
      Point x_j = (*it).node2().position();
      force_spring -= (*it).value().K * ((n.position() - x_j) / norm(n.position() - x_j)) * (norm(n.position() - x_j) - (*it).value().L);
    }
    return force_spring;
  }
};

/** Force function object that implements the damping force for HW2 #4.*/
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; // Quiet compiler warning
    return -n.value().vel * damping * Point(1);
  }
};

/** Force function object that returns 0 force (used as a placeholder) for HW2 #4. */
struct ZeroForce {
  /** Return the zero force (used as placeholder) applying to @a n at time @a t.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; (void)n;// Quiet compiler warning
    return Point(0);
  }
};


/** Force function object that sums three forces 
@param[in] three forces to be added
@return total force
*/
template <typename F1 = ZeroForce, typename F2 = ZeroForce, typename F3 = ZeroForce>
struct make_combined_force {
  /** Private Constructor */
  make_combined_force(F1 f1 = ZeroForce(), F2 f2 = ZeroForce() , F3 f3 = ZeroForce())
  : f1(f1), f2(f2), f3(f3) {}

  /* Atributes of the class. */
  F1 f1; // Force 1
  F2 f2; // Force 2
  F3 f3; // Force 3

  /** Return the total force applying to @a n at time @a t by adding the effect of individual forces.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
};

/** Plane constraint function object for HW2 #5a. */
struct PlaneConstraint {
  /** Add a constraint for a plane with the following properties:
   *     Plane z = −0.75.
   *     A node violates this constraint if x i · (0, 0, 1) < −0.75
   *   To fix a Node that violates this constraint:
   *     Set the position to the nearest point on the plane.
   *     Set the z-component of the Node velocity to 0
   */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t; // Quiet compiler warning
    // For every node check constraints.
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // Check constraint.
      if (dot(n.position(), Point(0, 0, 1)) < -0.75) 
      {
        n.position().z = - 0.75;
        n.value().vel *= Point(1, 1, 0);
      }
    }
  }
};

/** Sphere constraint function object for HW2 #5b. */
struct SphereConstraint {
  /* Add a constraint for a sphere with the following properties:
   *   Sphere center c = (0.5, 0.5, −0.5), sphere radius r = 0.15.
   *   A node violates this constraint if |x_i − c| < r.
   * To fix a Node that violates this constraint:
   *   Set the position to the nearest point on the surface of the sphere.
   *   Set the component of the velocity that is normal to the sphere’s surface to 0
   *   v_i = v_i − (v_i·R_i)R_i
   */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t; // Quiet compiler warning
    // For every node check for constraints.
    Point R;
    const double r = 0.15;
    const Point c = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // Check constraint.
      if (norm(n.position() - c) < r) 
      {
        n.position() = c + r / norm(n.position() - c) * (n.position() - c);
        R= (n.position() - c) / norm(n.position() - c);
        n.value().vel -= (dot(n.value().vel, R)) * R;
      }
    }
  }
};

/** No constraint function object used as placeholder for the constraints. */
struct NoConstraint  {
  /** Add a placeholder constraint that signifies no constraint.*/
  template <typename G>
  void operator()(G& g, double t) {
    (void)g; (void)t; // Quiet compiler warning
  }
};

/** Constraint function object that checks for all constraints for HW2 #5. */
template <typename C1 = NoConstraint, typename C2 = NoConstraint, typename C3 = NoConstraint>
struct constraints {
  /** Private Constructor */
  constraints(C1 c1 = NoConstraint(), C2 c2 = NoConstraint() , C3 c3 = NoConstraint())
  : c1(c1), c2(c2), c3(c3) {}

  /* Atributes of the class. */
  C1 c1; // Constraint1
  C2 c2; // Constraint2
  C3 c3; // Constraint3

  /** Return the total constraints applying to @a n at time @a t.*/
  template <typename G>
  void operator()(G& g, double t) {
    (void)t; // Quiet compiler warning
    // Check constraint.
    c1(g, t);
    c2(g, t);
    c3(g, t);
  }
};

/** Tear constraint function object for HW2 #6. */
struct TearConstraint {
  /* Add a constraint for a sphere with the following properties:
   *   Sphere center c = (0.5, 0.5, −0.5), sphere radius r = 0.15.
   *   A node violates this constraint if |x_i − c| < r.
   * To fix a Node that violates this constraint:
   *    Delete the node and all of its edges.
   */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t; // Quiet compiler warning
    // For every node check for constraints.
    const double r = 0.15;
    const Point c = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // Check constraint.
      if (norm(n.position() - c) < r) {
        g.remove_node(n);
      }
    }
  }
};
// iterates over nodes in the neighborhood
struct hoodnodes{
  const Point& cen_;
  Node n_;
  double rad_;
  hoodnodes(const Point& cen, Node n, double rad):
  cen_(cen), n_(n), rad_(rad) {}
  void operator()(Node n2)
  {
    Point r = cen_ -n2.position();
    double  l2 = normSq(r);
    if (n2!=n_ && l2<rad_)
    {
      n_.value().vel-= (dot(r, n_.value ().vel) / l2) * r; //remove r velocity component
    }
  }
};

struct nodecheck{
  SpaceSearcher<Node> search_;
  nodecheck(SpaceSearcher<Node> search): search_(search){}
  void operator()(Node n)
  {
    Edge e;
    const Point& center=n.position();
    double  radius2 = std::numeric_limits<double >::max();
    for(auto eit=n.edge_begin();eit!=n.edge_end();++eit)
    {
      e=*eit;
      radius2 = std::min(radius2 , normSq(e.node2 (). position () - center ));
    }
    radius2*=0.9;
    Box3D bb(Point(center + sqrt(radius2)), Point(Point(center - sqrt(radius2))));
    thrust::for_each(search_.begin(bb), search_.end(bb), hoodnodes(center, n, radius2));
  }
};
struct ParallelSelfCollisionConstraint{
  // template <typename G>
  void operator()(GraphType& g, double junk)
  {
    Box3D bigbox(Point(-3,-3,-3),Point(3,3,3));
    auto node2point=[](const Node& n){return n.position();};
    SpaceSearcher<Node> search(bigbox,g.node_begin(),g.node_end(),node2point);
    thrust::for_each(thrust::omp::par,g.node_begin(),g.node_end(),nodecheck(search));
  }
};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  for (auto n : nodes) {
    n.value().mass = 1.0/graph.num_nodes()/graph.num_nodes();
  }
  damping = 1.0 / float(nodes.size());
  Edge e = *(graph.edge_begin());
  L=e.length();
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    // Edge e = *it;
    (*it).value().L = (*it).length();
    (*it).value().K = 100.0/graph.num_nodes();
  }


  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
    double dt = 1.0/graph.num_nodes();
    double t_start = 0;
    double t_end = 5.0;

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
      // symp_euler_step(graph, t, dt, Problem1Force());
      /* symp_euler_step(graph, t, dt, Problem1Force(), NoConstraint()); */
      // symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce, DampingForce>(GravityForce(), MassSpringForce(), DampingForce()), NoConstraint());
      // symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce>(GravityForce(), MassSpringForce()), constraints<PlaneConstraint, SphereConstraint>(PlaneConstraint(), SphereConstraint()));
      symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce>
        (GravityForce(), MassSpringForce()), constraints<PlaneConstraint, SphereConstraint, ParallelSelfCollisionConstraint>
        (PlaneConstraint(), SphereConstraint(),ParallelSelfCollisionConstraint()));

        // Update viewer with nodes' new positions
      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
      viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
      if (graph.size() < 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
