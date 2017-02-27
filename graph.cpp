
// this is the number of generators of PSL2(Z)
// for groups with more generators, code should need only a few changes
#define LAB_CT 2


Graph::~Graph() {
   delete [] edges;
   delete [] edgeIniPN;
   delete [] edgeTerPN;
   delete [] vertIniHT;
   delete [] vertTerHT;
}


void Graph::Initialize(u32 Rep, u32 Rnu) {
   ep=Rep; //edge count
   nu=Rnu; //vertex count
   ds.Set(nu);
   edges = new GraphEdge [ep];
   edgeIniPN = new PrevNext [ep];
   edgeTerPN = new PrevNext [ep];
   vertIniHT = new HeadTail [nu];
   vertTerHT = new HeadTail [nu];
   checkPN = new PrevNext [nu];
   for (u32 k=0; k<ep; k++) {
      edgeIniPN[k].prev=-1; edgeIniPN[k].next=-1;
      edgeTerPN[k].prev=-1; edgeTerPN[k].next=-1;
      edges[k].ini=-1; edges[k].lab=-1; edges[k].ter=-1;
   }
   checkHT.head=-1; checkHT.tail=-1;
   for (u32 k=0; k<nu; k++) {
      checkPN[k].next=-1;   checkPN[k].prev=-1;
      vertIniHT[k].head=-1; vertIniHT[k].tail=-1;
      vertTerHT[k].head=-1; vertTerHT[k].tail=-1;
   }
}

Graph::Graph(u32 ep, u32 nu) {
   this->Initialize(ep,nu);
}

Graph::Graph(u32 mu, u32*S, u32*O) {
   this->Initialize(mu+mu,mu);
   u32 k=0;
   for(u32 n=0; n<mu; n++) {
      this->AddEdge(k,n,0,S[n]); k++;
      this->AddEdge(k,n,1,O[n]); k++;
   }
   assert(this->FoldedQ());
}



// return mu (and S and O)
// S and O are allocated
//   unless the Graph does not give a group of finite index
u32 Graph::ToGroup(u32*&Sout, u32*&Oout) {
   u32 k=0, j, h, t, f;
   u32* map=new u32 [nu];
   u32* P;
   for (j=0; j<nu; j++) {map[j]=-1;}
   for (j=0; j<nu; j++) {
      f=ds.Find(j);
      if (map[f]==-1) {map[f]=k;k++;}
      map[j]=map[f];
   }
   u32*S=new u32 [k];
   u32*O=new u32 [k];
   for (j=0; j<nu; j++) {
      if (ds.RootQ(j)) {
         h=vertIniHT[j].head; t=vertIniHT[j].tail;
         assert(h!=-1);assert(t!=-1);
         if (h==t) {goto InfiniteIndex;}
         assert(edgeIniPN[h].next==t);assert(edgeIniPN[t].prev==h);
         P = (edges[h].lab==0) ? S : O;
         P[map[ds.Find(edges[h].ini)]]=map[ds.Find(edges[h].ter)];
         P = (edges[t].lab==0) ? S : O;
         P[map[ds.Find(edges[t].ini)]]=map[ds.Find(edges[t].ter)];
      }
   }
   delete [] map;
   Sout=S;
   Oout=O;
   return k;

InfiniteIndex:
   delete [] map;
   delete [] S;
   delete [] O;
   Sout=Oout=NULL;
   return 0;

}


void Graph::Print() {
   printf("to check: ");
   PrintNodes(checkHT,checkPN);
   printf("\n");      
   printf("edges: ");
   for (u32 e=0; e<ep; e++) {
      if (edges[e].lab<LAB_CT) {
         assert(edges[e].ini!=-1); assert(edges[e].ter!=-1);
         printf("%d > %d > %d, ",edges[e].ini,/*'a'+*/edges[e].lab,edges[e].ter);
      }
   }
   printf("\n");
   for (u32 v=0; v<nu; v++) {
      if (ds.RootQ(v)) {
         printf("vert %d  ini: ",v);
         PrintNodes(vertIniHT[v],edgeIniPN);
         printf("  ter: ");
         PrintNodes(vertTerHT[v],edgeTerPN);
         printf("\n");
      }
   }
}

// Adds the edge k. Runs in O(LAB_CT) time.
void Graph::AddEdge(u32 k, u32 ini, u32 lab, u32 ter) {
   assert(edges[k].ini==-1 && edges[k].lab==-1 && edges[k].ter==-1);
   assert(ini<nu); assert(lab<LAB_CT); assert(ter<nu);
   edges[k].ini=ini; edges[k].lab=lab; edges[k].ter=ter;
   if (NodeFreeQ(ini,checkHT,checkPN)) {
      for (u32 e=vertIniHT[ini].head; e!=-1; e=edgeIniPN[e].next) {
         if (lab==edges[e].lab) {AppendNode(ini,checkHT,checkPN);}
      }
   }
   if (NodeFreeQ(ter,checkHT,checkPN)) {
      for (u32 e=vertTerHT[ter].head; e!=-1; e=edgeTerPN[e].next) {
         if (lab==edges[e].lab) {AppendNode(ter,checkHT,checkPN);}
      }
   }
   AppendNode(k,vertIniHT[ini],edgeIniPN);
   AppendNode(k,vertTerHT[ter],edgeTerPN);
}

// merge the two vertices v1 and v2. Runs in O(1) time.
// v1 and v2 don't have to be roots
// v2's edges are merged into v1's edges
void Graph::MergeVertices(u32 v1, u32 v2) {
   assert(v1!=-1); assert(v2!=-1);
   v1=ds.Find(v1); v2=ds.Find(v2);
   u32 v = ds.Merge(v1,v2);
   if (v!=v1) {
      JoinNodes(v,v1,vertIniHT,edgeIniPN);
      JoinNodes(v,v1,vertTerHT,edgeTerPN);
      if (NodeFreeQ(v,checkHT,checkPN)) {AppendNode(v,checkHT,checkPN);}
   }
   if (v!=v2) {
      JoinNodes(v,v2,vertIniHT,edgeIniPN);
      JoinNodes(v,v2,vertTerHT,edgeTerPN);
      if (NodeFreeQ(v,checkHT,checkPN)) {AppendNode(v,checkHT,checkPN);}
   }

}

// Check if the graph is folded.
// This could be replaced with a more robust check.
bool Graph::FoldedQ() {return checkHT.tail==-1;}

// Folds the graph. Runs in O(LAB_CT*ep*LogStar[ep]) time.
void Graph::Fold() {
   u32 e1, e2, e, v, u, w, l, i;
   u32 tab[LAB_CT];
   while ((v=checkHT.tail)!=-1) {
      //printf("loop start\n");
      //this->Print();getchar();
      DeleteNode(v,checkHT,checkPN);
      if (!ds.RootQ(v)) {continue;}
      for (i=0; i<LAB_CT; i++) {tab[i]=-1;}
      for (e=vertIniHT[v].head; e!=-1; e=edgeIniPN[e].next) {
         edges[e].ini=ds.Find(edges[e].ini); edges[e].ter=ds.Find(edges[e].ter);
         l=edges[e].lab; assert(l!=-1);
         if (tab[l]==-1) {tab[l]=e;}
         else {e1=tab[l]; e2=e; goto FoundIni;}
      }
      for (i=0; i<LAB_CT; i++) {tab[i]=-1;}
      for (e=vertTerHT[v].head; e!=-1; e=edgeTerPN[e].next) {
         edges[e].ini=ds.Find(edges[e].ini); edges[e].ter=ds.Find(edges[e].ter);
         l=edges[e].lab; assert(l!=-1);
         if (tab[l]==-1) {tab[l]=e;}
         else {e1=tab[l]; e2=e; goto FoundTer;}
      }
FoundNone:
      continue;
FoundIni:
      assert(v==edges[e1].ini); u=edges[e1].ter; w=edges[e2].ter;
      goto RemoveEdge;
FoundTer:
      assert(v==edges[e1].ter); u=edges[e1].ini; w=edges[e2].ini;
RemoveEdge:
      //printf("e1: %d  e2: %d\n",e1,e2);
      assert(ds.RootQ(v)); assert(ds.RootQ(u)); assert(ds.RootQ(w));
      DeleteNode(e2,vertIniHT[edges[e2].ini],edgeIniPN); edges[e2].ini=-1; 
      DeleteNode(e2,vertTerHT[edges[e2].ter],edgeTerPN); edges[e2].ter=-1; edges[e2].lab=-1;
      if (u!=w) {
         if (w==ds.Merge(u,w)) {swap(u,w);} //ensure u is the root
         JoinNodes(u,w,vertIniHT,edgeIniPN);
         JoinNodes(u,w,vertTerHT,edgeTerPN);
         if (NodeFreeQ(u,checkHT,checkPN)) {AppendNode(u,checkHT,checkPN);}
      }
      v=ds.Find(v);
      if (NodeFreeQ(v,checkHT,checkPN)) {AppendNode(v,checkHT,checkPN);}
   }
}




