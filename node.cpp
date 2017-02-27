void PrintNodes(HeadTail &ht, PrevNext*np) {
   printf("{");
   for (u32 h=ht.head; h!=-1; h=np[h].next) {
      printf("%d",h);
      if (np[h].next!=-1) {printf(", ");}
   }
   printf("}%d",ht.tail);
}

bool NodeElemQ(u32 c, HeadTail &ht, PrevNext*np) {
   return np[c].next!=-1 || ht.tail==c;
}
bool NodeFreeQ(u32 c, HeadTail &ht, PrevNext*np) {
   return np[c].next==-1 && ht.tail!=c;
}
void AppendNode(u32 c, HeadTail &ht, PrevNext*np) {
   u32 h=ht.head, t=ht.tail;
   if (h==-1) {ht.head=c; ht.tail=c; np[c].next=-1; np[c].prev=-1;}
   else {ht.tail=c; np[t].next=c; np[c].prev=t; np[c].next=-1;}
}
void DeleteNode(u32 c, HeadTail &ht, PrevNext*np) {
   u32 h=ht.head, t=ht.tail;
   u32 p=np[c].prev, n=np[c].next;
   if      (h==c&&c==t) {ht.head=-1;   ht.tail=-1;}
	else if (h==c&&c!=t) {ht.head=n;    np[n].prev=-1;}
	else if (h!=c&&c==t) {ht.tail=p;    np[p].next=-1;}
	else                 {np[p].next=n; np[n].prev=p;}
   np[c].next=-1; np[c].prev=-1;
}
void JoinNodes(u32 l1, u32 l2, HeadTail*ht, PrevNext*np) {
   u32 h1=ht[l1].head, t1=ht[l1].tail;
   u32 h2=ht[l2].head, t2=ht[l2].tail;
   if      (h2==-1) {return;}
   else if (h1==-1) {ht[l1].head=h2; ht[l1].tail=t2; ht[l2].head=-1;  ht[l2].tail=-1;}
   else  {np[t1].next=h2; np[h2].prev=t1; ht[l1].tail=t2; ht[l2].head=-1; ht[l2].tail=-1;}
}
