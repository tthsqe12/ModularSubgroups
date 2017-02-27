class Graph {
   u32 ep; //edge count
   u32 nu; //vertex count
   DisjointSet ds;
   GraphEdge*edges;
   PrevNext*edgeIniPN; // edgeIniPN[e].next = next edge that comes from initial vertex of e
   PrevNext*edgeTerPN; // dido for terminals
   HeadTail*vertIniHT; // vertexIniHT[v].head = head of the list of edges that come from vertex v
   HeadTail*vertTerHT; // dido for terminals
   PrevNext*checkPN;  // checkPN[v].next = next vertex after v in the list of vertices that need to be checked
   HeadTail checkHT;  // checkHT.head = head of list of list of vertices that need to be checked
   public:
   friend bool NodeElemQ(u32 c, HeadTail &ht, PrevNext*np);
   friend bool NodeFreeQ(u32 c, HeadTail &ht, PrevNext*np);
   friend void AppendNode(u32 c, HeadTail &ht, PrevNext*np);
   friend void DeleteNode(u32 c, HeadTail &ht, PrevNext*np);
   friend void JoinNodes(u32 l1, u32 l2, HeadTail*ht, PrevNext*np);
   void Initialize(u32 Rep, u32 Rnu);
   void Print();
   void AddEdge(u32 k, u32 ini, u32 lab, u32 ter);
   void MergeVertices(u32 v1, u32 v2);
   bool FoldedQ();
   void Fold();
   u32 ToGroup(u32*&Sout, u32*&Oout);
   Graph(u32 ep, u32 nu);
   Graph(u32 mu, u32*S, u32*O);
   ~Graph();
};