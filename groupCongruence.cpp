// test if a subgroup is congruence

bool Group::CongruenceQ(){
   assert(mu!=0);
   u32 *rel=new u32[24*mu];
   u32 k=gpCongruenceTest(mu,S,O,rel);
   delete [] rel;
   return k==0;
}

// return is index of the congruence closure
// Sout and Oout are allocated and filled in
u32 gpCongruenceClosure(u32 mu, u32*S, u32*O, u32*&Sout, u32*&Oout) {
   u32 *rel=new u32[24*mu];
   Graph graph(mu,S,O);
   u32 ct=gpCongruenceTest(mu,S,O,rel);
   for (u32 k=0; k<ct; k++) {
      for (u32 j=0; j<mu; j++) {
         graph.MergeVertices(j,*(rel+k*mu+j));
      }
   }
   graph.Fold();
   return graph.ToGroup(Sout,Oout);
}


// rel should have space for at least 24*n u32's (24 permutations)
// the return value is the number of relations that are not satisified
//  and the corresponding permutations are listed in rel, rel+n, ... .
// for example, if 3 is returned:
//    rel+0*n: some non-identity permutation
//    rel+1*n: some non-identity permutation
//    rel+2*n: some non-identity permutation
// The group is congruence iff zero is returned
u32 gpCongruenceTest(u32 n, u32*S, u32*O, u32*rel) {
   u32 ct=0;
   u32*cur=rel;
   u32*L,*R,*t1,*t2,*t3,*t4,*t5,*t6,*s,*a,*b,*r,*l,*abia,*lril;
      L=rel+8*n; R=rel+9*n;
      t1=rel+10*n; t2=rel+11*n; t3=rel+12*n; t4=rel+13*n; t5=rel+14*n; t6=rel+15*n;
      s=rel+16*n;  a=rel+17*n; b=rel+18*n; r=rel+19*n; l=rel+20*n; 
      abia=rel+21*n; lril=rel+22*n;
   pmMul(n,L,O,S);
   pmMul(n,R,O,L);
   u32 N=pmOrder(n,L,t1);
   u32 p=0, e=1, m=N;
      while ((m&1)==0) {p++; e=e+e; m=m/2;}
   i32 array[4] = {4,2,1,3};
   i32 inv5 = (1+(i32)e*(array[p%4]))/5;
   i32 inv2 = (1-(i32)m)/2;
   bool id;
   if (e==1) {
      pmMul(n,t1,R,R);
      pmPow(n,t2,L,-inv2,t4,t5);
      pmMul(n,t3,t1,t2);
      pmMul(n,cur,t3,t3,t3);
      if (!pmIdQ(n,cur)) {ct++;cur+=n;}
   }
   else if (m==1) {
      pmPow(n,t1,L,20,t4,t5);
      pmPow(n,t2,R,inv5,t4,t5);
      pmMul(n,t3,t1,t2);
      pmMul(n,t1,R,L,L,L,L);
      pmMulInv(n,s,t3,t1);
      pmMul(n,cur,s,S,s,S);
      if (!pmIdQ(n,cur)) {ct++;cur+=n;}

      pmMul(n,t5,R,R,R,R,R);
      pmMul(n,t1,t5,t5,t5,t5,t5);
      pmMul(n,t2,s,t1);
      pmMul(n,t3,R,s);
      id=pmMulInvIdQ(n,cur,t3,t2);
      if (!id) {ct++;cur+=n;}

      pmMul(n,t2,s,t5);
      pmMul(n,t3,t2,S);
      pmMul(n,cur,t3,t3,t3);
      if (!pmIdQ(n,cur)) {ct++;cur+=n;}
   }
   else {
		u32 c=e; while (c%m!=1){c+=e;} assert(c<N);
		u32 d=m; while (d%e!=1){d+=m;} assert(d<N);
      pmPow(n,a,L,c,t4,t5);
      pmPow(n,b,R,c,t4,t5);
      pmPow(n,l,L,d,t4,t5);
      pmPow(n,r,R,d,t4,t5);
      pmPow(n,t1,l,20,t4,t5);
      pmPow(n,t2,r,inv5,t4,t5);
      pmMul(n,t3,t1,t2);
      pmMul(n,t1,r,l,l,l,l);
      pmMulInv(n,s,t3,t1);

      pmMulInv(n,t1,a,b);
      pmMul(n,abia,t1,a);
      pmMulInv(n,t1,l,r);
      pmMul(n,lril,t1,l);

      pmMul(n,t1,a,r);
      pmMul(n,t2,r,a);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}

      pmMul(n,cur,abia,abia,abia,abia);
      if (!pmIdQ(n,cur)) {ct++;cur+=n;}

      pmMul(n,t1,abia,abia);
      pmInv(n,t2,b);
      pmMul(n,t3,t2,a);
      pmMul(n,t2,t3,t3,t3);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}

      pmPow(n,t2,a,-inv2,t4,t5);
      pmMul(n,t3,b,b,t2);
      pmMul(n,t2,t3,t3,t3);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}

      pmMul(n,t1,s,lril);
      pmMulInv(n,t2,lril,s);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}

      pmMul(n,t5,r,r,r,r,r);
      pmMul(n,t1,r,s);
      pmMul(n,t3,t5,t5,t5,t5,t5);
      pmMul(n,t2,s,t3);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}

      pmMul(n,t1,lril,lril);
      pmMul(n,t3,s,t5,lril);
      pmMul(n,t2,t3,t3,t3);
      id=pmMulInvIdQ(n,cur,t1,t2);
      if (!id) {ct++;cur+=n;}      
   }
   return ct;
}



