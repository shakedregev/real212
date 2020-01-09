/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <queue> 
#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;
using IncidentIter  = typename GraphType::incident_iterator;
/** Compares the distance of point x to point p with that of point y to p
* @param p The point which we are comparing distance to
* @param x,y References to nodes to see which is closer
* @ return true is x is close to p than y, false otherwise
**/
class distance{
public:
  Point p;
  distance(const Point& pp)
  : p(pp){};
  bool operator()(const NodeType x, const NodeType y)
  {
    return (normSq((x.position()-p))<(normSq(y.position()-p)));
  }
};
/** Compares the distance of point x to point p with that of point y to p
* @param longest longest path from a given root in the graph
* @param x References to a node to color
* @ return color for the node
* @ pre: x.value<=longest
**/
class heatmap{
public:
  int longest_;
  heatmap(int longest)
  : longest_(longest){};
  CME212::Color operator()(NodeType& x)
  {
    if (x.value()==-1)
      return CME212::Color::make_heat(0.0);
    return CME212::Color::make_heat(1.0-(float)x.value()/(float)(longest_));
    // return CME212::Color::make_heat((float)1);
  }
};

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */


NodeIter nearest_node(const GraphType& g, const Point& point)
{
  // HW1 #3: YOUR CODE HERE
  distance dist(point);
  return std::min_element(g.node_begin(),g.node_end(),dist);
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType& g, NodeType& root)
{
  // HW1 #3: YOUR CODE HERE
  std::queue <NodeType> n_queue;
  std::set <unsigned> n_visited;
  NodeIter niter=g.node_begin();
  while (niter!=g.node_end()) //initialize everything to -1
  {
    (*niter).value()=-1;
    ++niter;
  }
  int path_len=0;
  root.value()=path_len;
  n_queue.push(root);
  n_visited.insert(root.index());
  NodeType& curr=root;
  IncidentIter ince=root.edge_begin();
  while(!(n_queue.empty()))
  {
    //std::cout<<n_visited.size()<< " size" <<std::endl;
    curr=n_queue.front();
    n_queue.pop();
    path_len=curr.value();
    ince=curr.edge_begin();
    // int count=0;
    // std::cout<<n_visited.size()<< curr.degree() << " current"<<std::endl;
    while(ince!=curr.edge_end())
    {
      if(n_visited.find((*ince).node2().index())==n_visited.end())
      {
        (*ince).node2().value()=path_len+1;
        //std::cout<<"hi"<<std::endl;
        n_queue.push((*ince).node2());
        n_visited.insert((*ince).node2().index());
      }
      ++ince;
      // count++;
    }
    // std::cout<<n_visited.size()<< count<< " count" <<std::endl;
  }
  return curr.value();
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
      std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SFML_Viewer
      CME212::SFML_Viewer viewer;

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer

  // Center the view and enter the event loop for interactivity
      const GraphType &refgraph=graph;
      NodeType root=*(nearest_node(refgraph,Point(-1,0,1)));
      int longpath=shortest_path_lengths(graph,root);
      heatmap hm(longpath);
      auto nodemap=viewer.empty_node_map(graph);
      viewer.add_nodes(graph.node_begin(),graph.node_end(),hm,nodemap);
      viewer.add_edges(graph.edge_begin (), graph.edge_end(), nodemap );
      viewer.center_view();
      viewer.event_loop();

      return 0;
    }
