struct DisjointSetNode {
   u32 parent;
   u32 rank;
};

class DisjointSet {
   u32 size;
   DisjointSetNode* nodes;
   public:
   DisjointSet() {
      size=0; nodes=NULL;
   }
   DisjointSet(u32 n) {
      size=n;
      nodes = new DisjointSetNode[n];
      for (u32 i=0; i<n; i++) {
         nodes[i].rank=0; nodes[i].parent=i;
      }
   };
   void Set(u32 n) {
      if (nodes) {delete [] nodes;}
      size=n;
      nodes = new DisjointSetNode[n];
      for (u32 i=0; i<n; i++) {
         nodes[i].rank=0; nodes[i].parent=i;
      }
   }
   ~DisjointSet() {
      delete [] nodes;
      size=0;
   };
   void Print() {
      printf("{");
      for(u32 i=0;i<size;i++){
         printf("%d",nodes[i].parent);
         if(i+1<size) {printf(", ");}
      }
      printf("}");
   }
   u32 Find(u32 n) {
      if (nodes[n].parent==n) {return n;}
      else {return nodes[n].parent=Find(nodes[n].parent);}
   };
   u32 QFind(u32 n) {  // does not modify the forest
      if (nodes[n].parent==n) {return n;}
      else {return Find(nodes[n].parent);}
   };
   bool RootQ(u32 n) {return nodes[n].parent==n;};
   u32 Merge(u32 X, u32 Y) {
      u32 x=Find(X), y=Find(Y);
      if (x==y) {return -1;}
      if (nodes[x].rank>nodes[y].rank) {return nodes[y].parent=x;}
      else if (nodes[y].rank>nodes[x].rank) {return nodes[x].parent=y;}
      else {nodes[x].rank=nodes[x].rank+1; return nodes[y].parent=x;}
   }
};
