struct PrevNext {
   u32 prev;
   u32 next;   
};

struct HeadTail {
   u32 head;
   u32 tail;
};

struct GraphEdge {
   u32 ini;
   u32 ter;
   u32 lab;  //lab=-1 means no label, otherwise lab should be <LAB_CT
   u32 extra;
   GraphEdge() {;};
   GraphEdge(u32 ini1, u32 lab1, u32 ter1) {ini=ini1;lab=lab1;ter=ter1;};
};

void PrintNodes(HeadTail &ht, PrevNext*np);
bool NodeElemQ(u32 c, HeadTail &ht, PrevNext*np);
bool NodeFreeQ(u32 c, HeadTail &ht, PrevNext*np);
void AppendNode(u32 c, HeadTail &ht, PrevNext*np);
void DeleteNode(u32 c, HeadTail &ht, PrevNext*np);
void JoinNodes(u32 l1, u32 l2, HeadTail*ht, PrevNext*np);
