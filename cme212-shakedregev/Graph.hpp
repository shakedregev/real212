#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// Hw4 additions: thrust iterators
#include <thrust/system/omp/execution_policy.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/for_each.h>


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
	using node_value_type = V;

	using edge_value_type = E;

	using graph_type = Graph;

  /** Predeclaration of Node type. */
	class Node;
  /** Synonym for Node (following STL conventions). */
	using node_type = Node;

  /** Predeclaration of Edge type. */
	class Edge;
  /** Synonym for Edge (following STL conventions). */
	using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
	class NodeIterator;
  /** Synonym for NodeIterator */
	using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
	class EdgeIterator;
  /** Synonym for EdgeIterator */
	using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
	class IncidentIterator;
  /** Synonym for IncidentIterator */
	using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
	using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
	Graph() 
    // HW0: YOUR CODE HERE
	: nodes_gr(), adjacency(),ecount(0), ind2uid() {
// graph constructor: vector of nodes and edges
	}

  /** Default destructor */
	~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
	class Node : private totally_ordered<Node> {
	public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    /** Return this node's position. */
		const Point& position() const {
      // HW0: YOUR CODE HERE
			return mygraph->nodes_gr[nid].mypoint;
			//access point of relevant node
		}

		Point& position(){
      // modifiable position
			return mygraph->nodes_gr[nid].mypoint;
			//access point of relevant node
		}

    /** Return this node's index, a number in the range [0, graph_size). */
		size_type index() const {
      // HW0: YOUR CODE HERE
      // return index of relevant node
			return mygraph->nodes_gr[nid].nindex;
		}

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
/* Finds the degree of the node
@return degree
*/
		size_type  degree() const
		{
			return mygraph->adjacency[nid].size();
		}
/* Finds the first node connected to this node
@return the incident iterator initialized to this node
*/
		incident_iterator  edge_begin() const
		{
			return  incident_iterator(mygraph,nid, 0);
		}
/* Go in memory one past the last node connected to this node
@return the incident iterator initialized to this node
*/
		incident_iterator  edge_end() const
		{
			return  incident_iterator(mygraph,nid, this->degree());
		}
/* Finds the value of a node (we use it as distance from root node)
@return said value;
*/
		node_value_type& value()
		{
			return mygraph->nodes_gr[nid].value;
		}
/* Finds the value of a node (we use it as distance from root node)
@return said value as a constant;
*/
		const  node_value_type& value() const
		{
			// return const_cast<node_value_type&> (mygraph->nodes_gr[nid].value);
			return mygraph->nodes_gr[nid].value;
		}
		Node() {
      // HW0: YOUR CODE HERE
      //this is fine as is
		}


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
		bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // graph and node id must be equal
			return ((this->mygraph == n.mygraph)  && (this->nid == n.nid));
		}

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
		bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Graph must be less than or graph is equal and node id is less than
			if (this->mygraph<n.mygraph)
				return true;
			if ((this->mygraph==n.mygraph) && (this->nid <n.nid))
				return true;
			return false;
		}

	private:
    // Allow Graph to access Node's private member data and functions.
		friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
		Graph* mygraph;
		size_type nid;
		Node(const Graph* graph_, size_type nid_) 
		: mygraph(const_cast<Graph*>(graph_)), nid(nid_){
		} 
// initialize the graph and the id value of the node
	};

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
	size_type size() const {
    // HW0: YOUR CODE HERE
		return ind2uid.size();
// const time operation on a vector
	}

  /** Synonym for size(). */
	size_type num_nodes() const {
		return size();
	}

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
	Node add_node(const Point& position, const node_value_type& nodeval= node_value_type()) {
    // HW0: YOUR CODE HERE
    // create new node and give it index equal to the current size
		//works because this is 0 based
		grnodes new_node; 
		std::vector<gredges> blanklist;
		new_node.mypoint = position;
		new_node.value=nodeval;
		new_node.nindex = ind2uid.size();
		nodes_gr.push_back(new_node);
		adjacency.push_back(blanklist);
		// ind2uid.push_back(new_node.nindex);
		ind2uid.push_back(nodes_gr.size()-1);
		return Node(this, nodes_gr.size()-1);       
	}

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
	bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // check if the graph is the same, and large enough to have this id
		return (n.mygraph == this && nodes_gr[n.nid].nindex < num_nodes());
	}

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
	Node node(size_type i) const {
    // HW0: YOUR CODE HERE
		return Node(this, ind2uid[i]);

	}

  /**Remove specified node and incident edges
   * @param[in] n    Node to be removed.
   * @pre has_node(@a n) == True
   * @post has_node(@a n) == False
   * @post NodeIterators and EdgeIterators for @a n are invalidated.
   * @post IncidentIterators going through @a n's incident edges are invalidated.
   * @post If (e.node1() == @a n || e.node2() == @a n) -> e is removed.
   * @post new size() == old size() - 1
   * @return size_type  Original index of the node removed.
   * Complexity: O(num_nodes())
  */
	size_type remove_node(const Node& n) {
		if (!has_node(n)) {
			return 0;
		}

    // Delete edges adjacent to node n
		while (!adjacency[n.nid].empty()) {
			size_type dnode = adjacency[n.nid][0].eid;
			remove_edge(n, Node(n.mygraph, dnode));
		}

    // Deleting the element from index to nid vector
		auto iter = ind2uid.begin();
		while (iter != ind2uid.end()) {
			if (*iter == n.nid) {
				*iter = ind2uid.back();
				ind2uid.pop_back();
				break;
			} else {++iter;}
		}

    // decrementing the indices after the deletion
		for (unsigned i = nodes_gr[n.nid].nindex; i != ind2uid.size(); i++) {
			grnodes& in = nodes_gr.at(ind2uid[i]);
			in.nindex = i;
		}
		return *iter;
	}

  /**Remove specified node iterator
   * @param[in] n_it    NodeIterator to be removed.
   * @pre has_node(@a n) == True
   * @post has_node(@a n) == False
   * @post NodeIterators and EdgeIterators for @a n are invalidated.
   * @post IncidentIterators going through @a n's incident edges are invalidated.
   * @post If (e.node1() == @a n || e.node2() == @a n) then e is removed.
   * @post new size() == old size() - 1
   * @return node_iterator pointing to where original node was before removal.
   * Complexity: O(num_nodes())
  */

	node_iterator remove_node(node_iterator n_iter) {
		Node n = (*n_iter);
		size_type junk = remove_node(n);
		(void) junk;
		return n_iter;
		// return remove_node(n);
	}

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
	class Edge : private totally_ordered<Edge> {
	public:
    /** Construct an invalid Edge. */
		Edge() {
      // HW0: YOUR CODE HERE
      // this is fine 
		}

    /** Return a node of this Edge */
		Node node1() const {
      // HW0: YOUR CODE HERE
			return Node(this->mygraph, n1_id);      
		}

    /** Return the other node of this Edge */
		Node node2() const {
      // HW0: YOUR CODE HERE
			return Node(this->mygraph, n2_id);     
		}
		edge_value_type& value()
		{
			Edge ed;
			if (n2_id < n1_id) {
				ed = Edge(mygraph, n2_id, n1_id);
			} else {
				ed = Edge(mygraph, n1_id, n2_id);
			}
			size_type n = ed.node2().nid;

			for (auto it = mygraph->adjacency[n].begin();
				it != mygraph->adjacency[n].end(); ++it) {
				if ((*it).eid == ed.node1().nid)
					return (*it).value;
			}
			return mygraph->adjacency[ed.n1_id][ed.n2_id].value;
		}

		const  edge_value_type& value() const
		{
			Edge ed;
			if (n2_id < n1_id) {
				ed = Edge(mygraph, n2_id, n1_id);
			} else {
				ed = Edge(mygraph, n1_id, n2_id);
			}
			size_type n = ed.node2().nid;

			for (auto it = mygraph->adjacency[n].begin();
				it != mygraph->adjacency[n].end(); ++it) {
				if ((*it).eid == ed.node1().nid)
					return (*it).value;
			}
			return mygraph->adjacency[ed.n1_id][ed.n2_id].value;
		}


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
		bool operator==(const Edge& e) const {
      // (void) e;           // Quiet compiler warning
      // graph must be equal and both of the nodes in this must equal a node in e
			return ((this->mygraph==e.mygraph) && (((e.n1_id == n1_id) && (e.n2_id == n2_id)) || 
				((e.n1_id == n2_id) && (e.n2_id == n1_id))));
		}

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
		bool operator<(const Edge& e) const {
//if the graphs are not the same, ordering is determined by graph
  	//else ordering is determined by the smaller of the two indices
  	// if those are tied, the larger of the indices determines
			if (this->mygraph<e.mygraph)
				return true;
			if (this->mygraph>e.mygraph)
				return false;
			if ((e.n1_id < n1_id) && (e.n1_id<n2_id))
				return true;
			if ((e.n2_id < n1_id) && (e.n2_id<n2_id))
				return true;
			if ((e.n1_id == n1_id) && (e.n2_id<n2_id))
				return true;
			if ((e.n1_id < n1_id) && (e.n2_id==n2_id))
				return true;
			if ((e.n2_id == n1_id) && (e.n1_id<n2_id))
				return true;
			if ((e.n2_id < n1_id) && (e.n1_id==n2_id))
				return true;
			return false;
		}

		double length() const {
			return norm(node1().position() - node2().position());
		}

	private:
    // Allow Graph to access Edge's private member data and functions.
		friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
		Graph* mygraph;
		size_type n1_id; 
		size_type n2_id;
		Edge(const Graph* graph_, size_type n1_id_, size_type n2_id_)
		: mygraph(const_cast<Graph*>(graph_)), n1_id(n1_id_), n2_id(n2_id_){};  
// constructing and initializing the edge
	};

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
	size_type num_edges() const {
    // HW0: YOUR CODE HERE
		return ecount;
	}

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

	Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
		edge_iterator it = edge_begin();
		for (size_type j = 0; j < i; j++)
			++it;
		return *it;      
	}

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
	bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //check each edge to see if both of its nodes are equal to a or b
		for (auto iter = adjacency[a.nid].begin(); iter < adjacency[a.nid].end(); iter++) {
			if ((*iter).eid== b.nid) {
				return true;
			}
		}
		return false;
	}

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
	Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    // add edge with relevant properites
		if (has_edge(a,b)){
			return Edge(this, a.nid, b.nid);
		}else{
			// std::pair <size_type, size_type> new_edge; 
			// new_edge.n1_id = a.index(); 
			// new_edge.n2_id = b.index(); 
			// new_edge.eid=edges_gr.size();
			// new_edge.value=val;
			// edges_gr.push_back(new_edge);
			adjacency[a.index()].push_back({b.index(), edge_value_type()});
			adjacency[b.index()].push_back({a.index(), edge_value_type()});
			ecount++;
			return Edge(this, a.index(), b.index());     
		}
	}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
	void clear() {
    // HW0: YOUR CODE HERE
		nodes_gr.clear();
		ecount=0;
		adjacency.clear();
		ind2uid.clear();
	}

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
	// class NodeIterator : private totally_ordered<NodeIterator> {
	// public:
 //    // These type definitions let us use STL's iterator_traits.
 //    using value_type        = Node;                     // Element type
 //    using pointer           = Node*;                    // Pointers to elements
 //    using reference         = Node&;                    // Reference to elements
 //    using difference_type   = std::ptrdiff_t;           // Signed difference
 //    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

//     /** Construct an invalid NodeIterator. */
//     NodeIterator() {
//     }

//     // HW1 #2: YOUR CODE HERE
//     // Supply definitions AND SPECIFICATIONS for:
//     /* A dereferncer for nodes
//     @ return current node in the node iterator
//     */
//     Node operator*() const
//     {
//     	//return Node(graph_, graph_->nodes_gr[graph_->ind2uid[indx_]].nindex);
//     	return graph_->node(indx_);
//     }
//     /* An incrementor for node iterators
//     @ pre NodeIterator is not at the end (meaningless in this case)
//     @ return incremented node iterator
//     */
//     NodeIterator& operator++()
//     {
//     	indx_++;
//     	return *this;
//     }
//       A comparator for node iterators
//     @ return Returns true if graph and index of iterator are the same, false o.w.
       
//     bool operator==(const NodeIterator& n) const
//     {
//     	return ((n.graph_==this->graph_) && (this->indx_==n.indx_));
//     }

// private:
// 	friend class Graph;
// 	Graph* graph_;
// 	size_type indx_; 
// 	NodeIterator(const Graph* graph,size_type indx)
// 	: graph_(const_cast<Graph*>(graph)), indx_(indx){};
// };

//   // HW1 #2: YOUR CODE HERE
//   // Supply definitions AND SPECIFICATIONS for:
// /* Finds first index of node iterator
//    @ return iterator at said index
// */
// node_iterator node_begin() const
// {
// 	//return node_iterator(this,ind2uid.begin());
// 	return node_iterator(this,0);
// }
// /* Finds one past last index of node iterator
//    @ return iterator said index
// */
// node_iterator node_end() const
// {
// 	return node_iterator(this,ind2uid.size());
// }

struct id2node{ //converts ids to nodes
	id2node(const graph_type* g) : g_(g) {}
	Node operator ()(size_type id)
	{return Node(g_,id);}
	const graph_type* g_;
};

  struct NodeIterator : thrust::transform_iterator<id2node, std::vector<size_type>::const_iterator, Node> {
    using trans_it = thrust::transform_iterator<id2node, std::vector<size_type>::const_iterator, Node>;
    NodeIterator(){} 
    //conversion constructor .
    NodeIterator(const trans_it& t_iter) : trans_it{t_iter} {} 
  private:
    friend class Graph;
    NodeIterator(std::vector<size_type>::const_iterator it, const graph_type* g) : trans_it{it, id2node(g)} {}
  };

/* Finds first index of node iterator
   @ return iterator at said index
*/
node_iterator node_begin() const
{
	return node_iterator(ind2uid.begin(),this);
	// return node_iterator(this,0);
}
/* Finds one past last index of node iterator
   @ return iterator said index
*/
node_iterator node_end() const
{
	return node_iterator(ind2uid.end(),this);
	// return node_iterator(this,ind2uid.size());
}



  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
class IncidentIterator  : private totally_ordered<IncidentIterator> {
public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {

    }
    
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /* A dereferncer for edges
    @ return current edge in the incident iterator
    */
    Edge operator*() const{
    	return Edge(graph_, id1_, graph_->adjacency[id1_][id2_].eid);
    }
        /* An incrementor for incident iterators
    @ pre IncidentIterator is not at the end (meaningless in this case)
    @ return incremented incident iterator
    */
    IncidentIterator& operator++()
    {
    	if(id2_<graph_->adjacency[id1_].size())
    		id2_++;
    	return *this;
    }
    /* A comparator for incident iterators
    @ return Returns true if graph and both indices of iterator are the same, false o.w.
    */ 
    bool operator==(const IncidentIterator& n) const
    {
    	return ((n.graph_==this->graph_) && (this->id1_==n.id1_) && (this->id2_==n.id2_));
    }

private:
	friend class Graph;
	Graph* graph_;
	size_type id1_,id2_;
	IncidentIterator(const Graph* graph,size_type id1, size_type id2)
	: graph_(const_cast<Graph*>(graph)), id1_(id1), id2_(id2){};

};


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
class EdgeIterator : private totally_ordered<EdgeIterator> {
public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /* Dereferencer for edges in edge iterator
    @return the edge corressponding to the current place in the iterator
    */
    Edge operator*() const
    {
    	return Edge(graph_, graph_->adjacency[n1_id][eindex_].eid, n1_id);
    }
        /* incrementor for edges in edge iterator
    @return the iterator with edge index incremented
    */
        EdgeIterator& operator++()
    {
    	if (eindex_<graph_->adjacency[n1_id].size())
    		eindex_++;
    	get_edge();
    	return *this;
    }
        // Function that gets correct edge
    void get_edge() {
    	long adj = graph_->adjacency.size();
    	if (n1_id >= adj) return;
    	if (graph_->adjacency[n1_id].size() <= eindex_) {
    		n1_id++;
    		eindex_ = 0;
    	}

    	while ((n1_id < adj) && (graph_->adjacency[n1_id].size() == 0)) {
    		n1_id++;
    	}

    	while ((n1_id < adj) && (graph_->adjacency[n1_id][eindex_].eid < n1_id)) {
    		eindex_++;
    		if (eindex_ >= graph_->adjacency[n1_id].size()) {
    			n1_id++;
    			eindex_ = 0;
    		}
    	}
    	while ((n1_id < adj) && (graph_->adjacency[n1_id].size() == 0)) {
    		n1_id++;
    	}
    }


            /* comparator for edges in edge iterator
    @return true if the graph and edge index are the same, false o.w.
    */

    bool operator==(const EdgeIterator& n) const
    {
    	return ((n.graph_==this->graph_) && (this->eindex_==n.eindex_) && (this->n1_id==n.n1_id));
    }

private:
	friend class Graph;
    // HW1 #5: YOUR CODE HERE
	Graph* graph_;
	size_type n1_id;
	size_type eindex_;
	// EdgeIterator(const Graph* graph,size_type eindex, size_type id1,size_type id2)
	// : graph_(const_cast<Graph*>(graph)), eindex_(eindex), id1_(id1), id2_(id2){};
	EdgeIterator(const Graph* graph,size_type n1_id_, size_type eindex)
	: graph_(const_cast<Graph*>(graph)), n1_id(n1_id_), eindex_(eindex) {};
};

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
/* Go in memory to the first edge
@return the edge iterator initialized to this edge
*/
edge_iterator edge_begin() const
{
	return edge_iterator(this, 0,0);
}
/* Go in memory one past the edge node in the graph
@return the edge iterator initialized to this edge
*/
edge_iterator edge_end() const
{
	return edge_iterator(this, adjacency.size(), 0);
}


  /**Remove an edge from endpoints
  * @param[in] Node @a n1, first node of the edge to be removed
  * @param[in] Node @a n2, second node of the edge to be removed
  * @pre (has_node(n1) && has_node(n2)) == True
  * @post If (has_edge(n1, n2) == True) -> new num_edge() == old num_edge() - 1
  * @post If (has_edge(n1, n2) == False) -> new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return size_type  Index of the second node removed from the @a edges_ adjacency list
  * Complexity: O(num_nodes() + num_edges())
  */
size_type remove_edge(const Node& n1, const Node& n2) {
	size_type rmv_ind = 0;
	if (n2.nid >= adjacency.size() || n1.nid >= adjacency.size()) {
		return rmv_ind;
	}

	for (unsigned int i = 0; i != adjacency[n1.nid].size(); i++) {
		if (adjacency[n1.nid][i].eid == n2.nid) {
			adjacency[n1.nid].erase(adjacency[n1.nid].begin() + i);
			break;
		}
	}

	for (unsigned int i = 0; i != adjacency[n2.nid].size(); i++) {
		if (adjacency[n2.nid][i].eid == n1.nid) {
			adjacency[n2.nid].erase(adjacency[n2.nid].begin() + i);
			ecount--;
			rmv_ind = i;
			break;
		}
	}
	return rmv_ind;
}


  /**Remove an edge
  * @param[in] Edge e, edge to be removed
  * @pre @a e is a valid edge
  * @post If (has_edge(@a e.node1(), @a e.node2()) == True) -> new num_edge() == old num_edge()-1
  * @post If (has_edge(@a e.node1(), @a e.node2()) == False) -> new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return size_type  Index of the second node removed from the @a edges_ adjacency list
  * Complexity: O(num_nodes() + num_edges())
  */

size_type remove_edge(const Edge& e) {
	return remove_edge(e.node1(), e.node2());
}

  /**Remove an edge
  * @param[in] EdgeIterator @a e_it to be removed
  * @pre @a e_it is a valid edge iterator
  * @post If (has_edge(n1, n2) == True) -> new num_edge() == old num_edge() - 1
  * @post If (has_edge(n1, n2) == False) -> new num_edge() == old num_edge()
  * @post EdgeIterators, IncidentIterators invalidated where corresponding to the nodes
  * @return EdgeIterator  Points to the next valid edge in the @a edges_ adjacency list.
  * Complexity: O(num_nodes() + num_edges())
  */

edge_iterator remove_edge(edge_iterator e_it) {
	// size_type ind = remove_edge(*e_it);
	// (void) ind;
	// return e_it;
	return remove_edge(*e_it);
}


private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

// nodes have a point and an index
	// edges have two node indices
	struct grnodes
	{
		Point mypoint;
		size_type nindex;
		node_value_type value;
	};
	struct gredges
	{
		size_type eid;
		edge_value_type value;
		// size_type n1_id; 
		// size_type n2_id; 
	};
	std::vector<grnodes> nodes_gr;
	//std::vector<gredges> edges_gr;
	std::vector<std::vector<gredges>> adjacency;
	size_type ecount;
	std::vector<size_type> ind2uid;
};

#endif // CME212_GRAPH_HPP
