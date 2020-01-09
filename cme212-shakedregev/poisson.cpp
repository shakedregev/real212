/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"


// HW3: YOUR CODE HERE
// Define a GraphGraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
using GraphType = Graph<char,char>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;
using EdgeType = typename GraphType::edge_type;

/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  // (void) g; (void) bb;   //< Quiet compiler
  for(auto it = g.node_begin(); it != g.node_end(); ++it){
    if(bb.contains((*it).position())){
      g.remove_node(it);
    }
  }
  // return;
}
double f(Point x)
{
  return 5*cos(norm_1(x));
}

double g(Point x)
{
  if (std::abs(norm_inf(x)-1)<1e-10) return 0; //actual equality with floats rarely works
  if((norm_inf(x - Point(0.6, 0.6, 0)) < 0.2) || (norm_inf(x - Point(-0.6, -0.6, 0)) < 0.2) ||
    (norm_inf(x - Point(0.6, -0.6, 0)) < 0.2) || (norm_inf(x - Point(-0.6, 0.6, 0)) < 0.2))
  {
    return -0.2;
  }
  if(Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(x)) return 1;
  return -1;
}

struct GraphSymmetricMatrix {
  /* Symmetric Matrix constructor */
  GraphSymmetricMatrix(GraphType& gr)
  : gr_(gr){} 

  /* * Helper function to perform multiplication. Allows for delayed
   * evaluation of results.
   * Assign::apply(a, b) resolves to an assignment operation such as
   a += b, a -= b, or a = b.
   * @pre @a size(v) == size(w) */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v , VectorOut& w, Assign) const{
   for (auto it = gr_.node_begin(); it != gr_.node_end(); ++it) {
    double y = 0;
    for (auto it_2 = gr_.node_begin(); it_2 != gr_.node_end(); ++it_2) {
      if ((*it).index() == (*it_2).index() && g((*it).position()) != -1) {
        y += v[(*it_2).index()];
      }else if ((*it).index() != (*it_2).index() && (g((*it).position()) != -1 || g((*it_2).position()) != -1)) {
            // y += 0*v[(*it_2).index()];
            ; //no need to perform an operation
          }
          else {
            y += L((*it), (*it_2))*v[(*it_2).index()];
          }
        }
        Assign::apply(w[(*it).index()],y);
      }
    }

  /* * Matvec forwards to MTL ’s lazy mat_cvec_multiplier operator */
  template<typename Vector>
    mtl::vec::mat_cvec_multiplier <GraphSymmetricMatrix, Vector>
    operator*(const Vector& v) const {
      return {*this, v};
    }


    double L(NodeType n1, NodeType n2) const
    {
      if (n1.index()==n2.index())
        return -1.0*n2.degree();
      if (gr_.has_edge(n1,n2))
        return 1.0;
      return 0.0;
    }

  // dimension getter for matrix
    size_t dimen() const
    {
      return gr_.size();
    }

  private:
    GraphType gr_;

  };
/** The  number  of  elements  in the  matrix. */
  inline std:: size_t  size(const GraphSymmetricMatrix& A)
  {
    return A.dimen()*A.dimen();
  }

/** The  number  of rows in the  matrix. */
  inline std:: size_t  num_rows(const GraphSymmetricMatrix& A)
  {
    return A.dimen();
  }
/** The  number  of columns in the  matrix. */
  inline std:: size_t  num_cols(const GraphSymmetricMatrix& A)
  {
    return A.dimen();
  }

/* * Traits that MTL uses to determine properties of our GraphSymmetricMatrix . */
  namespace mtl{
    namespace ashape{

    /* * Define GraphSymmetricMatrix to be a non-scalar type . */
    template<>
      struct ashape_aux <GraphSymmetricMatrix> {
        typedef nonscal type;
      };
  } // end namespace ashape

  /* * GraphSymmetricMatrix implements the Collection concept
   * with value_type and size_type */
  template<>
  struct Collection <GraphSymmetricMatrix>{
    typedef double value_type;
    typedef unsigned size_type;
  };
} // end namespace mtl

template <class Real, class NodeColor, class PosFn, class OStream = std::ostream>
class visual_iteration : public itl::cyclic_iteration<Real> 
{
  typedef itl::cyclic_iteration<Real> super;
  void print_viewer()
  {
    /* Add elements to the graph */
    auto node_map = viewer_.empty_node_map(g);
    viewer_.clear();
    node_map.clear();
    viewer_.add_nodes(g.node_begin(), g.node_end(), NodeColor(u_), PosFn(u_),  node_map);
    viewer_.add_edges(g.edge_begin(), g.edge_end(), node_map);

    // Center view.
    viewer_.center_view();
  }

public:

  /* Define the constructor */
  template <class Vector>
  visual_iteration(const Vector& r0, int max_it_, Real tol_, GraphType& graph, CME212::SFML_Viewer& viewer, mtl::dense_vector<double>& u, OStream& out = std::cout, Real abstol_ = Real(0), int cycle_ = 50)
  : super(r0, max_it_, tol_, abstol_, cycle_,out), g(graph), viewer_(viewer), u_(u) {}
  
  //changes upon convergence
  bool finished() { return super::finished(); }

  //changes upon convergence
  template <typename T>
  bool finished(const T& r) 
  {
    bool finish= super::finished(r);
    print_viewer();
    return finish;
  }

private: 
  /* Define private attributes. */
  GraphType& g; 
  CME212::SFML_Viewer& viewer_;
  mtl::dense_vector<double>& u_;
};

struct PosFn{
  // constructor for position functor
  PosFn (mtl::dense_vector<double> u)
  : u_(u) {}
  Point operator() (const NodeType& n)
  {
    return Point(n.position().x,n.position().y,u_[n.index()]);
  }
  mtl::dense_vector<double> u_;
};

struct NodeColor{
  // constructor for color functor
  NodeColor (mtl::dense_vector<double> u)
  : u_(u) {
    min=*(std::min_element(u_.begin(),u_.end()));
    max=*(std::max_element(u_.begin(),u_.end()));

  }
  // operator that takes in a node and generates a color
  CME212::Color operator() (const NodeType& n)
  {
    if(min==max)
      return CME212::Color::make_heat(0.5);    
    return CME212::Color::make_heat((max-u_[n.index()])/(max-min)*(max-u_[n.index()])/(max-min));
  }
  mtl::dense_vector<double> u_;
  double min,max;
};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphGraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.
  mtl::dense_vector<double> u(graph.size(), 0.0); // initial guess for u
  mtl::dense_vector<double> b(graph.size(), 0.0); // initialize b
  GraphSymmetricMatrix A(graph); // Construct GraphSymmetricMatrix A 
  double sum_g;
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it){
    if(g((*it).position()) != -1){
      b[(*it).index()] = g((*it).position()); 
    }else{
      sum_g = 0;  
      for(auto edge_it = (*it).edge_begin(); edge_it != (*it).edge_end(); ++edge_it){ 
        if(g((*edge_it).node2().position()) != -1){
          sum_g += g((*edge_it).node2().position()); 
        }
      }
      b[(*it).index()] = h*h*f((*it).position()) - sum_g;
    }
  }
  CME212::SFML_Viewer viewer;
  // itl::cyclic_iteration<double> stop(b, 1000, 1.e-10, 0.0, 50); //this is for part 4
  visual_iteration<double, NodeColor, PosFn> stop(b, 1000, 1.e-10, graph, viewer, u);
// solve Au=b
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]()
  {
    itl::cg(A,u,b,stop);
  });

  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
