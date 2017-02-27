#include<stdio.h>
#include<vector>
#include<list>
#include<assert.h>

   //#include <iostream>
   
   typedef long long int64;
   typedef unsigned long long uint64;
   typedef int int32;
   typedef unsigned int uint32;

using namespace std;


uint32 gcd(uint32 a, uint32 b) { for(;;){ if(a==0){return b;} b%=a; if(b==0){return a;} a%=b;}}
uint64 gcd(uint64 a, uint64 b) { for(;;){ if(a==0){return b;} b%=a; if(b==0){return a;} a%=b;}}
uint32 lcm(uint32 a, uint32 b) { uint32 t=gcd(a,b); return t ? ((a/t)*b) : 0; }
uint64 lcm(uint64 a, uint64 b) { uint64 t=gcd(a,b); return t ? ((a/t)*b) : 0; }
int64 Sign(int64 x) { return (x != 0) || (x >> 63);}
int64 FlipSign(int64 x, int64 y) { return (x ^ (y >> 63)) - (y >> 63);}
int64 Abs(int64 x) { return FlipSign(x,x);}
void pmPrintCycles(uint32, uint32*);


class Matrix {
   int64 a,b,c,d;
   public:
   int64 get11()const {return a;};
   int64 get12()const {return b;};
   int64 get21()const {return c;};
   int64 get22()const{return d;};      
   Matrix (int64 aa, int64 bb, int64 cc, int64 dd){ a=aa; b=bb; c=cc; d=dd; };
   Matrix (){ a=0; b=0; c=0; d=0; };
   Matrix (Matrix & m){a=m.get11();b=m.get12();c=m.get21();d=m.get22();}
   Matrix ( const Matrix&  other):a(other.get11()){b=other.get12();c=other.get21();d=other.get22();}
   void Set(int64 aa, int64 bb, int64 cc, int64 dd){ a=aa; b=bb; c=cc; d=dd; };
   int64 CuspRank() {return 8*Abs(a+c)+Abs(b+d);};
   void Print() {printf("{{ %d , %d } , { %d , %d }}",a,b,c,d);};
   void RMulBy(Matrix m) { // this = this.m
      int64 aa = a*m.get11() + b*m.get21();
      int64 bb = a*m.get12() + b*m.get22();
      int64 cc = c*m.get11() + d*m.get21();
      int64 dd = c*m.get12() + d*m.get22();
      a = aa; b = bb; c = cc; d = dd;
   };
   void LMulByInverse(Matrix m) { // this = inverse(m) . this
      int64 aa = a*m.get22() - c*m.get12();
      int64 bb = b*m.get22() - d*m.get12();
      int64 cc = c*m.get11() - a*m.get21();
      int64 dd = d*m.get11() - b*m.get21();
      a = aa; b = bb; c = cc; d = dd;
   };
   void Inverse() {
      swap(a,d); b = -b; c = -c;
   };
};


void Print(vector<Matrix> v){
   for(uint32 i=0; i<v.size(); i++) {
      printf("%d: ",i); v[i].Print(); printf(" \n");
   }
}
bool pOkQ(uint32 n, uint32*p) {
   
   for(uint32 i=0;i<n;i++)
   { if( p[i]>=n) return false;
     else { 
      for(uint32 j=0;j<n;j++)
      if(p[j]==p[i]&&(i!=j)) return false;
      }
   }
   return true;
};
   
class Group {
	uint32 mu;
	uint32* S;
	uint32* O;
	public:
	bool MemberQ( Matrix m);
   Group( const Group &other){
      mu=other.mu;
   	S = new uint32 [mu]; for (uint32 i=0;i<mu;i++){S[i]=other.S[i];};
		O = new uint32 [mu]; for (uint32 i=0;i<mu;i++){O[i]=other.O[i];};
	};
   Group(uint32 muu, uint32* SS, uint32* OO) {
	   mu = muu;
		S = new uint32 [muu]; for (uint32 i=0;i<muu;i++){S[i]=SS[i];};
		O = new uint32 [muu]; for (uint32 i=0;i<muu;i++){O[i]=OO[i];};
	};
   Group(bool (* f) (Matrix));
   //Group(bool (* f) (uint32 , Matrix),uint32 );
   vector<Matrix> Cosets();
   //vector<Matrix> Generators();
   bool SameQ(Group &that);
   bool ConjugateQ(Group & that);
   void Print();
   friend bool pOkQ(uint32 , uint32 *);
   Group Conjugate( uint32 p []){
      assert(pOkQ(mu,p));
      Group cj(*this);
      for(uint32 i=0;i<mu;i++)
      {cj.O[p[i]]=p[O[i]];cj.S[p[i]]=p[S[i]];};
      return cj;
   };
   void GPrint(vector<Matrix> v){
   for(uint32 i=0; i<v.size(); i++) {
      printf("%d: ",i); v[i].Print(); printf(" \n");
      }
   };
   bool HSameQ(Group &other,uint32 start);
   bool CongruenceQ();

	~Group(){//printf("calling group destructor\n");
      delete[] S; delete[] O;};
};

struct CosetInd {
   Matrix mRep{1,0,0,1};
   uint32 ind=-1;
   void Print() {mRep.Print();printf("  %d",ind);};
};





Group::Group(bool (* f) (Matrix)){
   if ((*f)(Matrix(0,-1,1,0)) && (*f)(Matrix(1,1,0,1))) {
      mu=1; S = new uint32 [1] {0}; O = new uint32 [1] {0}; return;
   }
   if ((*f)(Matrix(0,-1,1,1)) && (*f)(Matrix(1,-1,0,1))) {
      mu=2; S = new uint32 [2] {1,0}; O = new uint32 [2] {0,1}; return;
   }
   uint32 p, q, k=3;
   vector <uint32>s; s.push_back(-1); s.push_back(-1); s.push_back(-1);
   vector <uint32>o; o.push_back(1); o.push_back(2); o.push_back(0);
   list <CosetInd>l;
   CosetInd c;
   c.mRep.Set(1,0,0,1); c.ind=0; l.push_back(c);
   c.mRep.Set(1,-1,1,0); c.ind=1; l.push_back(c);
   c.mRep.Set(0,-1,1,-1); c.ind=2; l.push_back(c);
   Matrix m(0,0,0,0);

loop:
   for(list<CosetInd>::iterator it = l.begin(); it != l.end();) {      
      m.Set(0,-1,1,0); m.RMulBy((*it).mRep); m.LMulByInverse((*it).mRep);
      if ((*f)(m)) {
         p = (*it).ind; assert(s[p]==-1); s[p] = p;
         it++; l.erase(prev(it)); continue;
      }
      m.Set(0,1,-1,-1); m.RMulBy((*it).mRep); m.LMulByInverse((*it).mRep);
      if ((*f)(m)) {
         p = (*it).ind; assert(s[p]==-1); s[p] = k; s.push_back(p); o.push_back(k); k += 1;
         it++; l.erase(prev(it)); continue;
      }
      for(list<CosetInd>::iterator jt = l.begin(); jt != it; jt++){
         m.Set(0,-1,1,0); m.RMulBy((*jt).mRep); m.LMulByInverse((*it).mRep);
         if((*f)(m)) {
            p = (*it).ind; q = (*jt).ind;
            assert(s[p]==-1); assert(s[q]==-1); s[p] = q; s[q] = p;
            l.erase(jt); it++; l.erase(prev(it)); it--; break;          
         }
      }
      it++;
   }

   int64 d, bd = 1000000000000;
   list<CosetInd>::iterator bi = l.begin();
   if (l.size() == 0) {goto done;}
   if (k>32000) {goto failed;}
   for(list<CosetInd>::iterator it = l.begin(); it != l.end(); it++) {      
      d = (*it).mRep.CuspRank(); if (d < bd) {bd = d; bi = it;}
   }
   p = (*bi).ind;
   o.push_back(k+1); o.push_back(k+2); o.push_back(k);
   assert(s[p]==-1); s[p] = k; s.push_back(p); s.push_back(-1); s.push_back(-1);
   c.mRep.Set(1,1,0,1); c.mRep.RMulBy((*bi).mRep); c.ind=k+1; l.push_back(c);
   c.mRep.Set(1,0,1,1); c.mRep.RMulBy((*bi).mRep); c.ind=k+2; l.push_back(c);
   l.erase(bi);
   k += 3;
   goto loop;

done:
   mu = k;
   S = new uint32 [k];
   O = new uint32 [k];
   assert(k==s.size()); assert(k==o.size());
   for (uint32 i = 0; i<k; i++) {S[i]=s[i];};
   for (uint32 i = 0; i<k; i++) {O[i]=o[i];};
   return;

failed:
   mu = 1;
   S = new uint32 [1] {0};
   O = new uint32 [1] {0};
   return;
};

bool Group::MemberQ(Matrix  m) {
   int64 m11 = m.get11(), m12 = m.get12(), m21 = m.get21(), m22 = m.get22();
   int64 n11, n12, n21, n22, a1, a2, b1, b2;
   uint32 k=0;
   for (uint32 it=0;it<1000;it++) {
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
   return false;
}


vector<Matrix> Group::Cosets(){
   vector <Matrix> cs;
   cs.resize(mu);
   cs[0].Set(1,0,0,1);
   if(mu == 1) { return cs;}
   if(mu == 2) {cs[1].Set(0,-1,1,0); return cs;}
   cs[1].Set(1,-1,1,0);
   cs[2].Set(0,-1,1,-1);

   uint32 p, q, k=3;
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
   
   int64 d, bd = 1000000000000;
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

bool Group::HSameQ(Group &other,uint32 start) {
   uint32 i1 = 0, i2 = start,ni1=0, ni2=0;
   if (mu != other.mu){ return false;}
   if (mu < 3) {return true;}

   bool ret = false;

   uint32*Ores1 = new uint32 [mu];
   uint32*Ores2 = new uint32 [mu];
   uint32*map = new uint32 [mu];
   for(uint32 i = 0; i< mu; i++) { Ores1[i] = Ores2[i] =map[i]= -1;}
   for(uint32 i=0;i<3*mu;i++){
      //printf("%d\n",i+1);
      if(map[i1]==-1) {map[i1]=i2;}
      else if(map[i1]!=i2) {goto done;}

      if(i1==O[i1]){ //printf("target o has valance 1\n");//target o has valence 1
         if(i2!=other.O[i2]) {goto done;}
         i1=S[i1];
         i2=other.S[i2];
      }
      else { //printf("target o has valance 3\n");//target o has valence 3
         ni1=O[i1];ni2=other.O[i2];
         //printf("need back up\n");// if need back up 
         if(Ores1[i1]!=-1&&Ores1[i1]!=i1) {
            if(Ores2[i2]==-1||Ores2[i2]==i2) {goto done;}
            else i1=S[i1];i2=other.S[i2];
         }
         else {//printf("no need back up\n");//no need back up 
            //printf("target e has valence 1\n");//target e has valence 1
            if(ni1==S[ni1]) { 
               if(ni2!=other.S[ni2]) {goto done;}
               else {
                  Ores1[i1]=Ores1[ni1]=Ores1[O[ni1]]=ni1;
                  Ores2[i2]=Ores2[ni2]=Ores2[other.O[ni2]]=ni2;
                  i1=ni1;i2=ni2;
               }
             }
            else {//printf("goto next\n");//go to next 
               if(map[ni1]==-1) map[ni1]=ni2;
               else if(map[ni1]!=ni2) {goto done;}
               Ores1[i1]=Ores1[ni1]=Ores1[O[ni1]]=ni1;
               Ores2[i2]=Ores2[ni2]=Ores2[other.O[ni2]]=ni2;
               i1=S[ni1];i2=other.S[ni2];
            }

         }
      }
      //printf("return to beginning\n");
      if(i1==0&&(Ores1[i1]==-1||Ores1[i1]==i1))
         if(!(i2==start&&(Ores2[i2]==-1||Ores2[i2]==i2))) {goto done;}
         else {ret=true; goto done;}
   }

done:
   //printf("done\n");
   delete[] Ores1; delete[] Ores2; delete [] map;
   //printf("delete done\n");
   return ret;
   

}

bool Group::SameQ(Group &other){ 
return this->HSameQ(other,0);

}

bool Group::ConjugateQ(Group &other){ 
   bool ret=false;
   for(uint32 i=0;i<mu;i++) 
   ret=ret||HSameQ(other,i);
   return ret;
}


/*bool Group::SameQ(Group &other) {
   uint32 i1 = 0, i2 = 0,ni1=0, ni2=0;
   if (mu != other.mu){ return false;}
   if (mu < 3) {return true;}

   bool ret = false;

   uint32*Ores1 = new uint32 [mu];
   uint32*Ores2 = new uint32 [mu];
   uint32*map = new uint32 [mu];
   for(uint32 i = 0; i< mu; i++) { Ores1[i] = Ores2[i] =map[i]= -1;}
   for(uint32 i=0;i<3*mu;i++){
      printf("%d\n",i+1);
      if(map[i1]==-1) {map[i1]=i2;}
      else if(map[i1]!=i2) {goto done;}

      if(i1==O[i1]){ printf("target o has valance 1\n");//target o has valence 1
         if(i2!=other.O[i2]) {goto done;}
         i1=S[i1];
         i2=other.S[i2];
      }
      else { printf("target o has valance 3\n");//target o has valence 3
         ni1=O[i1];ni2=other.O[i2];
         printf("need back up\n");// if need back up 
         if(Ores1[i1]!=-1&&Ores1[i1]!=i1) {
            if(Ores2[i2]==-1||Ores2[i2]==i2) {goto done;}
            else i1=S[i1];i2=other.S[i2];
         }
         else {printf("no need back up\n");//no need back up 
            printf("target e has valence 1\n");//target e has valence 1
            if(ni1==S[ni1]) { 
               if(ni2!=other.S[ni2]) {goto done;}
               else {
                  Ores1[i1]=Ores1[ni1]=Ores1[O[ni1]]=ni1;
                  Ores2[i2]=Ores2[ni2]=Ores2[other.O[ni2]]=ni2;
                  i1=ni1;i2=ni2;
               }
             }
            else {printf("goto next\n");//go to next 
               if(map[ni1]==-1) map[ni1]=ni2;
               else if(map[ni1]!=ni2) {goto done;}
               Ores1[i1]=Ores1[ni1]=Ores1[O[ni1]]=ni1;
               Ores2[i2]=Ores2[ni2]=Ores2[other.O[ni2]]=ni2;
               i1=S[ni1];i2=other.S[ni2];
            }

         }
      }
      printf("return to beginning\n");
      if(i1==0&&(Ores1[i1]==-1||Ores1[i1]==i1))
         if(!(i2==0&&(Ores2[i2]==-1||Ores2[i2]==i2))) {goto done;}
         else {ret=true; goto done;}
   }

done:
   printf("done\n");
   delete[] Ores1; printf("delete 1 done\n"); delete[] Ores2; delete [] map;
   printf("delete done\n");
   return ret;
   

}
*/

void pmPow (uint32 n,uint32*z, uint32*p, int32 pow){
     uint32* q=new uint32[n]; 
     uint32* q2=new uint32[n]; 
     for(uint32 i=0;i<n;i++) {q[i]=p[i];z[i]=i;}
     if(pow<0) 
     {pow=-pow; for(uint32 i=0;i<n;i++) q[p[i]]=i;}
     while(pow>0){ for(uint32 i=0;i<n;i++)q2[i]=q[i];
     if(pow%2==1) {for(uint32 i=0;i<n;i++) z[i]=q[z[i]];}
     pow=pow/2;
     for(uint32 i=0;i<n;i++) {q[i]=q2[q[i]];}
     }
     delete[] q; delete[] q2;
  

}

void pmMul(uint32 n, uint32*p, uint32*q, uint32*r) {
   for(uint32 i=0;i<n;i++){p[i]=q[r[i]];}
   return;
}


bool pmIdQ(uint32 n, uint32 *p){
   bool ret=false;
   for(uint32 i=0;i<n;i++) ret=ret||(p[i]==i);
   return ret;

}

uint32 pmOrder(uint32 n,uint32*p){
   uint32 order=1;
   uint32 *q=new uint32[n];
   uint32 *r=new uint32[n];
   pmPow(n,q,p,1);
   while(!pmIdQ(n,q)) {
   pmMul(n,r,p,q);
   for(uint32 i=0;i<n;i++) q[i]=r[i];
   order++;}
   return order; 
}
bool PermutationOK(uint32 n, uint32*p) {
   uint32 i, j, k;
   bool * visited = new bool[n];
   for(i=0;i<n;i++){visited[i]=false;};
   for(i=0;i<n;i++){
      if( (p[i] >= n) || visited[p[i]]){ delete[] visited; return false;}
      visited[p[i]]=true;
   };
   delete[] visited;
   return true;
}



uint32 pmFixedCount(uint32 n, uint32*p) {
   uint32 ct = 0;
   for(uint32 i=0;i<n;i++) {ct += (p[i]==i);}
   return ct;
};

uint32 pmCycleCount(uint32 n, uint32*p) {
   uint32 i, j, k, ct = 0;
   bool * visited = new bool[n];
   assert(PermutationOK(n,p));
   for(i=0;i<n;i++){visited[i]=false;};
   for(k=0;k<n;k++) {
      if (!visited[k]) {
         visited[k]=true;
         for(j=p[k]; j!=k; j=p[j]) {
            visited[j]=true;
         }
         ct++;
      }
   }
   delete[] visited;
   return ct;
}



void pmPrintCycles(uint32 n, uint32*p) {
   uint32 i, j, k;
   bool * visited = new bool[n];
   assert(PermutationOK(n,p));
   for(i=0;i<n;i++){visited[i]=false;};
   for(k=0;k<n;k++) {
      if (!visited[k]) {
         printf("(%d",k+1); visited[k]=true;
         for(j=p[k]; j!=k; j=p[j]) {
            printf(" %d",j+1); visited[j]=true;
         }
         printf(")");
      }
   }
   printf("\n");
   delete[] visited;
   return;
}
uint32 ModInv(uint32 n,uint32 a){
   assert((a!=0)&&(n>1)); 
   for(uint32 i=1;i<n;i++)
   if(i*a%n==1) return i;
}
bool Group::CongruenceQ(){
   uint32 * r=new uint32[mu];
   uint32 * v=new uint32[mu];
   uint32 * t=new uint32[mu];
   //v=s.o.o.s 
   pmPow(mu,v,O,2);
   pmMul(mu,r,v,S);
   pmMul(mu,v,S,r);
   //l=s.v^(-1) R=S.v^(-2)
   uint32 * L=new uint32[mu];
   uint32 * R=new uint32[mu];
   pmPow(mu,r,v,-1); //r=v^(-1)  
   pmMul(mu,L,S,r);
   pmPow(mu,r,v,-2);//r=v^(-2)
   pmMul(mu,R,S,r);
   uint32 N=pmOrder(mu,L);
   pmPow(mu,r,L,-ModInv(mu,2));
   pmPow(mu,v,R,2);
   pmMul(mu,t,v,r);
   pmPow(mu,v,t,3);
   return pmIdQ(mu,v);
}
void Group::Print(){
   uint32*t = new uint32[mu];
   pmMul(mu,t,O,S);
   printf("Group[(mu, e2, e3, ei) = (%d %d %d %d)\n    S = ",
            mu, 
            pmFixedCount(mu,S),
            pmFixedCount(mu,O),
            pmCycleCount(mu,t));
   pmPrintCycles(mu,S);
   printf("\n    O = ");
   pmPrintCycles(mu,O);
   printf("\n    T = ");
   pmPrintCycles(mu,t);
   printf("]");
   delete [] t;
}


bool GammaN(Matrix m) {
uint32 n = 3;
   int64 a = m.get11(), b = m.get12(), c = m.get21(), d = m.get22();
   return (c%n==0)&&(b%n==0)&&((((a-1)%n==0)&&((d-1)%n==0))||(((a+1)%n==0)&&((d+1)%n==0)));
}
bool Gamma1N(Matrix m) {
uint32 n = 15;
   int64 a = m.get11(), b = m.get12(), c = m.get21(), d = m.get22();
   return (c%n==0)&&((((a-1)%n==0)&&((d-1)%n==0))||(((a+1)%n==0)&&((d+1)%n==0)));
}
bool Gamma0N(Matrix m) {
uint32 n = 11;
   int64 a = m.get11(), b = m.get12(), c = m.get21(), d = m.get22();
   return (c%n==0);
}

 

//bool (gg)(Matrix) {return true;};

//bool gg(Matrix m) {return GammaN(2,m);}

int main()
{
   printf("Hello!\n");
   uint32 S[18] = {17, 1, 8, 4, 3, 9, 7, 6, 2, 5, 10, 12, 11, 14, 13, 16, 15, 0};
   uint32 O[18] = {17, 3, 1, 2, 6, 4, 5, 15, 13, 11, 9, 10, 8, 12, 7, 14, 0, 16};
   printf("\ng: \n");
   Group g(18, S, O);
   g.Print();
   printf("\n");


   Group h(Gamma0N);
   printf("h:\n");
   h.Print(); printf("\n"); 

   Matrix m(4,1,7,2);
   printf("%d\n",h.MemberQ(m));
   m.Set(2,1,701,351);
   printf("%d\n",h.MemberQ(m));

   vector<Matrix> cosets = h.Cosets();
   printf("h.Cosets():  {\n");
   Print(cosets);
   printf("\n");
   
   Group h2(GammaN);
   printf("h same h ?: %d\n",h.SameQ(h));
  
   getchar();
   h2.Print();
   getchar();
   uint32 p[12]={1,11,3,2,6,5,4,7,8,9,10,0};
   Group h3(h2.Conjugate(p));h3.Print();getchar();
   printf("? same: %d\n",h2.SameQ(h3));
   printf("? conjugate: %d\n",h2.ConjugateQ(h3));
   printf("permutation power\n");
   uint32* r=new uint32[18];pmPow(12,r,p,2);
   pmPrintCycles(12,r);
   printf("\n");
   for(uint32 i=0;i<4;i++){
   printf("permutation power%d\n",i);
   pmPow(18,r,O,i);
   pmPrintCycles(18,r);}
   printf("permutation composition\n");
   pmMul(18,r,O,O);
   pmPrintCycles(18,r);
   printf("%d\n",pmOrder(18,O)); 
   printf("Congruence test:\n");       
   printf("%d\n",h.CongruenceQ());    
      
   return 0;
   
};   
   
   
   
   
   