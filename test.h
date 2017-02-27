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


void PrintCycles(uint32*, uint32);

class Matrix {
   int64 a,b,c,d;
   public:
   int64 get11()const {return a;};
   int64 get12()const {return b;};
   int64 get21()const {return c;};
   int64 get22()const{return d;};      
   Matrix (int64 aa, int64 bb, int64 cc, int64 dd){ a=aa; b=bb; c=cc; d=dd; };
   Matrix (Matrix & m){a=m.get11();b=m.get12();c=m.get21();d=m.get22();}
   Matrix ( const Matrix&  other):a(other.get11()){b=other.get12();c=other.get21();d=other.get22();}
   void Set(int64 aa, int64 bb, int64 cc, int64 dd){ a=aa; b=bb; c=cc; d=dd; };
   void Print() {printf("{{ %d , %d } , { %d , %d }}",a,b,c,d);};
   RMulBy(Matrix m) { // this = this.m
      int64 aa = a*m.get11() + b*m.get21();
      int64 bb = a*m.get12() + b*m.get22();
      int64 cc = c*m.get11() + d*m.get21();
      int64 dd = c*m.get12() + d*m.get22();
      a = aa; b = bb; c = cc; d = dd;
   };
   LMulByInverse(Matrix m) { // this = inverse(m) . this
      int64 aa = a*m.get22() - c*m.get12();
      int64 bb = b*m.get22() - d*m.get12();
      int64 cc = c*m.get11() - a*m.get21();
      int64 dd = d*m.get11() - b*m.get21();
      a = aa; b = bb; c = cc; d = dd;
   };
   Inverse() {
      swap(a,d); b = -b; c = -c;
   };
};

   
class Group {
	uint32 mu;
	uint32* S;
	uint32* O;
	public:
	bool MemberQ( Matrix m);
	Group(uint32 muu, uint32* SS, uint32* OO) {
	   mu = muu;
		S = new uint32 [muu]; for (uint32 i=0;i<muu;i++){S[i]=SS[i];};
		O = new uint32 [muu]; for (uint32 i=0;i<muu;i++){O[i]=OO[i];};
	};
   Group(bool (* f) (Matrix));
   Group(bool (* f) (uint32 , Matrix),uint32 );
   void Print();
	~Group(){delete[] S; delete[] O;};
};

struct Coset {
   Matrix mRep{1,0,0,1};
   uint32 index=-1;
   void Print() {mRep.Print();printf("  %d",index);};
};



int64 sign(int64 x) {
   if (x>0) {return 1;}
   else if (x<0) {return -1;}
   else {return 0;}
};



Group::Group(bool (* f) (Matrix)){

   if ((*f)(Matrix(0,-1,1,0)) && (*f)(Matrix(1,1,0,1))) {
      mu=1;
      S = new uint32 [1] {0};
      O = new uint32 [1] {0};
      return;
   }
   if ((*f)(Matrix(0,-1,1,1)) && (*f)(Matrix(1,-1,0,1))) {
      mu=2;
      S = new uint32 [2] {1,0};  /* (01)*/
      O = new uint32 [2] {0,1};
      return;
   }
   uint32 k=3;
   vector <uint32>s; s.push_back(-1); s.push_back(-1); s.push_back(-1);
   vector <uint32>o; o.push_back(1); o.push_back(2); o.push_back(0);
   list <Coset>l;
   Coset c;
   c.mRep.Set(1,0,0,1); c.index=0; l.push_back(c);
   c.mRep.Set(1,-1,1,0); c.index=1; l.push_back(c);
   c.mRep.Set(0,-1,1,-1); c.index=2; l.push_back(c);
   uint32 

   l.front().Print();printf("\n");


loop:
   Matrix m;
   uint32 p, q;

   for(list<Coset>::iterator it = l.begin(); it < l.end(); it++) {      
      printf("Checking even pairings");
      m.Set(0,-1,1,0); m.RMulBy((*it).mRep); m.LMulByInverse((*it).mRep);
      if ((*f)(m)) {                                                       printf("even pairing");
         p = (*it).index; assert(s[p]==-1); s[p] = p;
         l.erase(it); it--; continue;
      }
      printf("Checking free pairings");
      m.Set(0,1,-1,-1); m.RMulBy((*it).mRep); m.LMulByInverse((*it).mRep);
      if ((*f)(m)) {                                                       printf("odd pairing");
         p = (*it).index; assert(s[p]==-1); s[k+1] = p; s.push_back(p); k += 1;
         l.erase(it); it--; continue;
      }
      printf("Checking free pairings");
      for(list<Coset>::iterator jt = l.begin(); jt != it; jt++){
         m.Set(0,-1,1,0); m.RMulBy((*jt).mRep); m.LMulByInverse((*it).mRep);
         if((*f)(m)) {                                                     printf("free pairing %p %p", it, jt);      
            p = (*it).index; q = (*jt).index;
            assert(s[p]==-1); assert(s[q]==-1); s[p] = q; s[q] = p;
            l.erase(jt); l.erase(it);  it -= 2; continue             
         }
      }
   }

   if (l.size() == 0) {goto done;}
   int64 d, bestd = 1000000000000;
   list<Coset>::iterator besti = l.begin();
   for(list<Coset>::iterator it = l.begin(); it < l.end(); it++) {      
      d = 24*abs((*it).mRep.get11()+(*it).mRep.get21())
         + 3*abs((*it).mRep.get21()+(*it).mRep.get22())
         - sign((*it).mRep.get21()+(*it).mRep.get22());
      if (d < bestd) {bestd = d; besti = it}
   }
   p = (*besti).index;
   o.push_back(k+2); o.push_back(k+3); o.push_back(k+1)
   assert(s[p]==-1); s[p] = k+1; s.push_back(p); s.push_back(-1); s.push_back(-1);
   m.Set(1,1,0,1); m.RMulBy((*besti).mRep); l.push_back(m);
   m.Set(1,0,1,1); m.RMulBy((*besti).mRep); l.push_back(m);
   l.erase(besti);
   k += 3;
   goto loop;

done:
   mu = k;
   S = new uint32 [k];
   O = new uint32 [k];
   assert(k==s.size());
   assert(k==o.size());
   for (uint32 i = 0; i<k; i++) {S[i]=s[i]; O[i]=o[i];}
   return;
};


void PermutationOK(uint32*p, uint32 n) {
   uint32 i, j, k;
   bool * visited = new bool[n];
   for(i=0;i<n;i++){visited[i]=false;};
   for(i=0;i<n;i++){
      if(!(p[i]<n ) || visited[p[i]]){return false;}
      visited[p[i]]=true;
   };
   return true;
};

void PrintCycles(uint32*p, uint32 n) {
   uint32 i, j, k;
   bool * visited = new bool[n];
   assert(PermutationOK(p,n));
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
}

void Group::Print(){
   printf("mu: %d\n",mu);
   printf("S: "); PrintCycles(S,mu);
   printf("\nO: "); PrintCycles(O,mu);
   printf("\n");
}


bool GammaN(Matrix m) {
uint32 n = 2;
   int64 a = m.get11();
   int64 b = m.get12();
   int64 c = m.get21();
   int64 d = m.get22();
   return (c%n==0)&&(b%n==0)&&((((a-1)%n==0)&&((d-1)%n==0))||(((a+1)%n==0)&&((d+1)%n==0)));
}

 
bool Group::MemberQ( Matrix  m) {
   int64 n11, m11 = m.get11();
   int64 n12, m12 = m.get12();
   int64 n21, m21 = m.get21();
   int64 n22, m22 = m.get22();
   int64 a1,a2,b1,b2;
   uint32 k=0;
   for (uint32 it=0;it<1000;it++) {
      n11=m11; n12=m12; n21=m21; n22=m22;
      //printf("{{ %d , %d } , { %d , %d }}\n",m11,m12,m21,m22);
      if (m21 == 0) {
         if (m12 == 0) { return k == 0; }
         else if ( (m12 ^ m22) >= 0) { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[k]; /*printf("1 O ");*/ }
         else { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("1 S ");*/ }
      }
      else if (m22 == 0) {
         if (m11 == 0) { return S[k] == 0; }
         else if ( (m21 ^ m11) >= 0) { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[k]; /*printf("2 O ");*/ }
         else { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("2 S ");*/ }
      }
      else {
         if (m21>0) {a1=m11;a2=m21;} else {a1=-m11;a2=-m21;}
         if (m22>0) {b1=m12;b2=m22;} else {b1=-m12;b2=-m22;}
         if (a1+b1<0) { m11 = n21; m12 = n22; m21 = -n11; m22 = -n12; k = S[k]; /*printf("3 S ");*/ }
         else if (a1+b1<a2+b2) { m11 = n11-n21; m12 = n12-n22; m21 = n11; m22 = n12; k = O[O[k]]; /*printf("3 OO");*/ }
         else { m11 = n21; m12 = n22; m21 = n21 - n11; m22 = n22 - n12; k = O[k]; /*printf("3 O ");*/ }
      }
   };
   printf("MemberQ max iterations");
   return false;
};

//bool (gg)(Matrix) {return true;};

//bool gg(Matrix m) {return GammaN(2,m);}

int main()
{
   printf("Hello!\n");
   uint32 S[18] = {17, 1, 8, 4, 3, 9, 7, 6, 2, 5, 10, 12, 11, 14, 13, 16, 15, 0};
   uint32 O[18] = {17, 3, 1, 2, 6, 4, 5, 15, 13, 11, 9, 10, 8, 12, 7, 14, 0, 16};
   PrintCycles(O,18);
   printf("\n g: \n");
   Group g(18, S, O);
   g.Print();
   printf("\n");

   bool (Group:: * gg) (Matrix)=&Group::MemberQ;
/* bool (*gt) (Matrix)=GammaN(5,);  
   Group h(( g.*gg));
*/

   Group h(GammaN);

   printf("h:\n");
   h.Print();

   printf("\n h print done \n");

   vector<Coset> mlist;
   Coset c;
   c.mRep.Set(2,4,5,6);
   c.index=1;
   printf("%lu\n",c.mRep.get11());
   mlist.push_back(c);
   printf("%lu\n",mlist[0].mRep.get11());
   c.mRep.Set(3,4,5,6);
   mlist.push_back(c);
   c.mRep.Set(4,4,5,6);
   mlist.push_back(c);
   printf("%lu\n",mlist[2].mRep.get11());

   
   Matrix m(6,1,17,3);
   if((g.*gg)(m)) {printf("true");} else {printf("false");};	  
   return 0;

};
   
   
   /*
   sint64 m11,m12,m21,m22
   
   m11=; m12=;
   m21=;
   */
   
   
   
   
   
   