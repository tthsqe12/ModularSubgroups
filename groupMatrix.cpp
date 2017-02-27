
// various functions for dealing with the matrices related to a group


// test if a matrix is in a group
bool Group::MemberQ(Matrix  m) {
   i64 m11 = m.get11(), m12 = m.get12(), m21 = m.get21(), m22 = m.get22();
   i64 n11, n12, n21, n22, a1, a2, b1, b2;
   u32 k=0;
   for (u32 it=0;it<1000;it++) {
      n11=m11; n12=m12; n21=m21; n22=m22;//      printf("{{ %d , %d } , { %d , %d }}\n",m11,m12,m21,m22);
      if (m21 == 0) {
         if (m12 == 0) { return k == 0; }
         else if ( (m12 ^ m22) >= 0) { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[O[k]]; /*printf("1 OO ");*/ }
         else { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("1 S ");*/ }
      }
      else if (m22 == 0) {
         if (m11 == 0) { return S[k] == 0; }
         else if ( (m21 ^ m11) >= 0) { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[O[k]]; /*printf("2 OO ");*/ }
         else { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("2 S ");*/ }
      }
      else {
         if (m21>0) {a1=m11;a2=m21;} else {a1=-m11;a2=-m21;}
         if (m22>0) {b1=m12;b2=m22;} else {b1=-m12;b2=-m22;}
         if (a1+b1<0) { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("3 S ");*/ }
         else if (a1+b1<a2+b2) { m11 = n11-n21; m12 = n12-n22; m21 = n11; m22 = n12; k = O[k]; /*printf("3 O");*/ }
         else { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[O[k]]; /*printf("3 OO ");*/ }
      }
   };
   printf("MemberQ max iterations");
   assert(false);
   return false;
}


struct CosetInd {
   Matrix mRep{1,0,0,1};
   u32 ind=-1;
   void Print() {mRep.Print();printf("  %d",ind);};
};



// given a function f for testing if a given matrix is in the group,
//  construct the representation of the group
// If the index gets too large or the function f is bad, the function fails
Group::Group(bool (* f) (Matrix)){
   if ((*f)(Matrix(0,-1,1,0)) && (*f)(Matrix(1,1,0,1))) {
      mu=1; S = new u32 [1] {0}; O = new u32 [1] {0}; return;
   }
   if ((*f)(Matrix(0,-1,1,1)) && (*f)(Matrix(1,-1,0,1))) {
      mu=2; S = new u32 [2] {1,0}; O = new u32 [2] {0,1}; return;
   }
   vector <u32>s;
   vector <u32>o;
   vector <u32>li;
   vector <Matrix>lm;
   u32 si, ei, i, j, p, q, k;
   if ((*f)(Matrix(1,-1,1,0))) {
      s.push_back(1);s.push_back(0);s.push_back(-1);s.push_back(-1);
      o.push_back(0);o.push_back(2);o.push_back(3);o.push_back(1);
      lm.push_back(Matrix(1,1,0,1)); li.push_back(2);
      lm.push_back(Matrix(1,0,1,1)); li.push_back(3);
      si=0; ei=1; k=4;
   }
   else {
      s.push_back(-1);s.push_back(-1);s.push_back(-1);
      o.push_back(1);o.push_back(2);o.push_back(0);
      lm.push_back(Matrix(1,0,0,1)); li.push_back(0);
      lm.push_back(Matrix(1,-1,1,0)); li.push_back(1);
      lm.push_back(Matrix(0,-1,1,-1)); li.push_back(2);
      si=0; ei=2; k=3;
   }
   Matrix m, mi;
Loop:

   printf("list:   si: %d  ei: %d \n",si+1,ei+1);
   for (i=0; i<=ei; i++) {
      printf("  %d ",li[i]+1);
      lm[i].Print();
      printf("\n");
   }

   for (i=si; i<ei+1; i++) {
      m.Set(0,-1,1,0); m.RMulBy(lm[i]); m.LMulByInverse(lm[i]);
      if ((*f)(m)) {
         p=li[i]; assert(s[p]==-1); s[p]=p;
         lm[i]=lm[ei];
         li[i]=li[ei];
         ei--; i--; continue;
      }
      m.Set(0,-1,1,1); m.RMulBy(lm[i]); m.LMulByInverse(lm[i]);
      if ((*f)(m)) {
         p=li[i]; assert(s[p]==-1); s[p]=k; s.push_back(p); o.push_back(k); k++;
         lm[i]=lm[ei];
         li[i]=li[ei];
         ei--; i--; continue;
      }
      for( j=0; j<i; j++){
         m.Set(0,-1,1,0); m.RMulBy(lm[j]); m.LMulByInverse(lm[i]);
         if((*f)(m)) {
            p=li[i]; q=li[j]; assert(s[p]==-1); assert(s[q]==-1); s[p]=q; s[q]=p;
            lm[j]=lm[i-1]; lm[i-1]=lm[ei]; lm[i]=lm[ei-1];
            li[j]=li[i-1]; li[i-1]=li[ei]; li[i]=li[ei-1];
				i-=2; ei-=2; break;          
         }
      }

   }
   i64 d, bd=999999999;
   if (ei==-1) {goto Done;}
   if (k>320) {goto Failed;}

	i=0; for (j=0; j<=ei; j++) {d=lm[j].CuspRank(); if (d<bd) {bd=d;i=j;}}
	mi=lm[i]; p=li[i];
   o.push_back(k+1); o.push_back(k+2); o.push_back(k);
   assert(s[p]==-1); s[p]=k; s.push_back(p); s.push_back(-1); s.push_back(-1);
	si=ei;
	lm[i]=lm[ei]; li[i]=li[ei];
   m.Set(1,1,0,1); m.RMulBy(mi); li[ei]=o[s[p]]; lm[ei]=m;
	ei++;
   m.Set(1,0,1,1); m.RMulBy(mi);
   if (ei+1>lm.size()) {
      assert(ei==lm.size()); assert(ei==li.size());
      li.push_back(o[o[s[p]]]); lm.push_back(m);
   }
   else {
      li[ei]=o[o[s[p]]]; lm[ei]=m;
   }
	k=k+3;
   goto Loop;
   
Done:
   mu = k; S = new u32 [k]; O = new u32 [k];
   assert(k==s.size()); assert(k==o.size());
   for (i=0; i<k; i++) {S[i]=s[i];O[i]=o[i];};
   return;

Failed:
   mu = 0; S = NULL; O = NULL;
   return;
}



// enumerate the cosets
vector<Matrix> Group::Cosets(){
  vector <Matrix> cs;
  cs.resize(mu);
  cs[0].Set(1,0,0,1);
  if(mu == 1) { return cs;}
  if(mu == 2) {cs[1].Set(0,-1,1,0); return cs;}
  cs[1].Set(1,-1,1,0);
  cs[2].Set(0,-1,1,-1);

  u32 p, q, k=3;
  list <CosetInd>l;
  CosetInd c;
  c.mRep.Set(1,0,0,1); c.ind=0; l.push_back(c);
  c.mRep.Set(1,-1,1,0); c.ind=1; l.push_back(c);
  c.mRep.Set(0,-1,1,-1); c.ind=2; l.push_back(c);
  Matrix m(0,0,0,0);

loop:
  for(list<CosetInd>::iterator it = l.begin(); it != l.end();) { 
    p = (*it).ind;      
    if (S[p] == p) {
      it++; l.erase(prev(it)); continue;
    }
    if (O[S[p]] == S[p]) {
      m.Set(0,-1,1,0); m.RMulBy((*it).mRep); cs[S[p]] = m;
      it++; l.erase(prev(it)); continue;
    }
    for(list<CosetInd>::iterator jt = l.begin(); jt != it; jt++){
      q = (*jt).ind;
      if(S[p] == q) {   assert(S[q] == p);
          l.erase(jt); it++; l.erase(prev(it)); it--; break;          
        }
    }
    it++;
  }

   //printf("coset:  \n"); GPrint(cs);
   //printf("l:  \n ");   
   //for(list<CosetInd>::iterator jt = l.begin(); jt != l.end(); jt++){ (*jt).Print(); printf(" \n");   }


   //getchar();
   
   i64 d, bd = 1000000000000;
   list<CosetInd>::iterator bi = l.begin();
   if (l.size() == 0) {goto done;}
   for(list<CosetInd>::iterator it = l.begin(); it != l.end(); it++) {      
      d = (*it).mRep.CuspRank(); if (d < bd) {bd = d; bi = it;}
   }
   p = (*bi).ind;
   c.mRep.Set(0,-1,1,0); c.mRep.RMulBy((*bi).mRep); cs[S[p]] = c.mRep;
   //printf("old: ");c.Print();printf(" \n");

   c.mRep.Set(1,1,0,1); c.mRep.RMulBy((*bi).mRep); cs[O[S[p]]] = c.mRep;
   c.ind = O[S[p]]; l.push_back(c);
   //printf("new: ");c.Print();printf(" \n");
   
   c.mRep.Set(1,0,1,1); c.mRep.RMulBy((*bi).mRep); cs[O[O[S[p]]]] = c.mRep;
   c.ind = O[O[S[p]]]; l.push_back(c);
   //printf("new: ");c.Print();printf(" \n");
   

   l.erase(bi);

   goto loop;

done:   
   return cs;
}




void gpCosetsGens(u32 mu, u32*S, u32*O, u32*cs, u32*gen) {


   return;

}



// construct the group from a list of generators
//  function fails if generators generate an infinite index subgroup
Group::Group(vector<Matrix> matl) {mu = gpFromGens(matl,S,O);}

// returns mu and allocates a new S and O (unless mu=0)
u32 gpFromGens(vector<Matrix> matl, u32*&S, u32*&O) {
   u32 nu=1;
   vector <GraphEdge> sm;
   for (u32 i=0;i<matl.size();i++){
      i64 m11=matl[i].get11(), m12=matl[i].get12(), m21=matl[i].get21(), m22=matl[i].get22();
      i64 n11, n12, n21, n22, a1, a2, b1, b2;
      u32 toAp=0, sta=-1, end=0, AP_NONE=0, AP_S=1, AP_O=2, AP_OO=3, AP_DONE=4; 
      while (toAp!=AP_DONE) {
         n11=m11; n12=m12; n21=m21; n22=m22;
         if (m21==0) {
            if (m12==0) {toAp=AP_DONE;}
            else if ((m12^m22)>=0) {m11=n21; m12=n22; m21=n21-n11; m22=n22-n12; toAp=AP_O;}
            else { m11=n21; m12=n22; m21=-n11; m22=-n12; toAp=AP_S;}
         }
         else if (m22==0) {
            if (m11==0) {m11=1; m12=0; m21=0; m22=1; toAp=AP_S;}
            else if ((m21^m11)>=0) {m11=n21; m12=n22; m21=n21-n11; m22=n22-n12; toAp=AP_O;}
            else {m11=n21; m12=n22; m21=-n11; m22=-n12; toAp=AP_S;}
         }
         else {
            if (m21>0) {a1=m11;a2=m21;} else {a1=-m11;a2=-m21;}
            if (m22>0) {b1=m12;b2=m22;} else {b1=-m12;b2=-m22;}
            if (a1+b1<0) {m11=n21; m12=n22; m21=-n11; m22=-n12; toAp=AP_S;}
            else if (a1+b1<a2+b2) {m11=n11-n21; m12=n12-n22; m21=n11; m22=n12; toAp=AP_OO;}
            else {m11=n21; m12=n22; m21=n21-n11; m22=n22-n12; toAp=AP_O;}
         }
         sta=end; if (m12==0&&m21==0) {end=0;} else {end=nu;nu++;}
         if (toAp==AP_S) {
            sm.push_back(GraphEdge(sta,0,end));
            sm.push_back(GraphEdge(end,0,sta));
         }
         else if (toAp==AP_OO) {
            sm.push_back(GraphEdge(sta,1,end));
            if (sta!=end) {
               sm.push_back(GraphEdge(end,1,nu));
               sm.push_back(GraphEdge(nu,1,sta));
               nu++;
            }
         }
         else if (toAp==AP_O) {
            sm.push_back(GraphEdge(end,1,sta));
            if (sta!=end) {
               sm.push_back(GraphEdge(nu,1,end));
               sm.push_back(GraphEdge(sta,1,nu));
               nu++;
            }
         }
      }
   }
   u32 ep = sm.size();
   if (ep==0) {/*printf("$InfiniteIndex\n");*/ S=O=NULL; return 0;}
   Graph graph(ep,nu);
   for (u32 i=0; i<ep; i++) {
      graph.AddEdge(i,sm[i].ini,sm[i].lab,sm[i].ter);
   }
//   graph.Print();
   graph.Fold();
   return graph.ToGroup(S,O);
};




