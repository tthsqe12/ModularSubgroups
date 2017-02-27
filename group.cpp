
Group::~Group() {if (mu!=0) {delete[] S; delete[] O;}};

Group Group::Conjugate(u32*p){
      assert(pmValidQ(mu,p));
      Group cj(*this);
      for(u32 i=0;i<mu;i++)
      {cj.O[p[i]]=p[O[i]];cj.S[p[i]]=p[S[i]];};
      return cj;
   };


// *************** graphs ******************

Group::Group(Graph &gr) {mu=gr.ToGroup(S,O);}







void Group::Standardize() {
   if (mu!=0) {
      u32*map = new u32 [mu];
      u32*Sn = new u32 [mu];
      u32*On = new u32 [mu];
      u32 k=gpStandardize(mu,S,O,map,Sn,On);
      assert(k!=0);
      delete [] map;
      delete [] S;
      delete [] O;
      S=Sn; O=On;
   }
}

u32 gpStandardize_f(u32 k, u32 n, u32*map, u32*S, u32*O) {
   if (map[n]==-1) {
      map[n]=k++;
      k=gpStandardize_f(k,O[n],map,S,O);
      k=gpStandardize_f(k,S[n],map,S,O);
   }
   return k;
}

u32 gpStandardize(u32 mu, u32*S, u32*O, u32*map, u32*Sn, u32*On) {
   assert(mu!=0);
   u32 i;
   i=0; do {map[i]=-1;} while (++i<mu);
   u32 k=gpStandardize_f(0,0,map,S,O);
   if (k<mu) {return 0;}
   assert(k==mu);
   i=0; do {Sn[map[i]]=map[S[i]];On[map[i]]=map[O[i]];} while (++i<k);
   return k;
}





// **************** printing ********************

void Group::Print(){
   if (mu==0) {printf("$BadGroup"); return;}
   u32*t = new u32[mu];
   u32*T = new u32[mu];
   pmMul(mu,T,O,S);
   printf("Group[mu = %d, e2 = %d, e3 = %d, ei = %d,\n   S = ",
            mu, 
            pmFixedCount(mu,S),
            pmFixedCount(mu,O),
            pmCycleCount(mu,T,t));
   pmPrintCycles(mu,S);
   printf(",\n   O = ");
   pmPrintCycles(mu,O);
   printf(",\n   T = ");
   pmPrintCycles(mu,T);
   printf("]");
   delete [] t;
   delete [] T;
}

