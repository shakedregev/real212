/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>>
{
public:
  // Get all of the iterator traits and make them our own
  using value_type        = typename std::iterator_traits<It>::value_type;
  using pointer           = typename std::iterator_traits<It>::pointer;
  using reference         = typename std::iterator_traits<It>::reference;
  using difference_type   = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
  : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
    while ((it_!=end_) && (p_(*it_) == false))
      ++it_;
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /* dereferencer for the filter iterator
  @return derferenced place in filter iterator
  */
  value_type operator*() const
  {
    return *it_;
  }
    /* incrementor for the filter iterator
  @return the new incremented filter iterator
  */
  filter_iterator& operator++()
  {
    ++it_;
    while ((p_(*it_) == false) && (it_!=end_))
      ++it_;
    return *this;
  }
  /* comparator for the filter iterators
  @return true if the current index in the iterators is the same, false o.w.
  */
  bool operator==(const filter_iterator& iter) const
  {
    return (iter.it_==it_);
  }


private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
 const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/*This predicate keeps only non-isolated nodes with even degree
@return true if condition is met, false otherwise
*/
struct Parity {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return ((n.degree() % 2 == 0) && (n.degree()!=0));
  }
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph<int>;
  using NodeType  = typename GraphType::node_type;

  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

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

  // HW1 #4: YOUR CODE HERE
  // Use the filter_iterator to plot an induced subgraph.
      auto nodemap=viewer.empty_node_map(graph);
      // auto it0=make_filtered(graph.node_begin(),graph.node_end(),SlicePredicate());
      // auto itend=make_filtered(graph.node_end(),graph.node_end(),SlicePredicate());
      auto it0=make_filtered(graph.node_begin(),graph.node_end(),Parity());
      auto itend=make_filtered(graph.node_end(),graph.node_end(),Parity());
      viewer.add_nodes(it0,itend,nodemap);
      viewer.add_edges(graph.edge_begin (), graph.edge_end(), nodemap );
  // Center the view and enter the event loop for interactivity
      viewer.center_view();
      viewer.event_loop();

      return 0;
    }
