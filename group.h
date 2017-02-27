class Group {
	u32 mu;
	u32* S;
	u32* O;
	public:
	bool MemberQ( Matrix m);
   Group( const Group &other){
      mu=other.mu;
   	S = new u32 [mu]; for (u32 i=0;i<mu;i++){S[i]=other.S[i];};
		O = new u32 [mu]; for (u32 i=0;i<mu;i++){O[i]=other.O[i];};
	};
   Group(u32 muu, u32* SS, u32* OO) {
	   mu = muu;
		S = new u32 [muu]; for (u32 i=0;i<muu;i++){S[i]=SS[i];};
		O = new u32 [muu]; for (u32 i=0;i<muu;i++){O[i]=OO[i];};
	};
   Group(Graph &gr);
   Group(vector<Matrix> matl);
   Group(bool (* f) (Matrix));
   void Standardize();
   vector<Matrix> Cosets();
   bool SameQ(Group &that);
   bool ConjugateQ(Group & that);
   void Print();
   friend bool pValidQ(u32 , u32 *);
   Group Conjugate(u32*p);
   void GPrint(vector<Matrix> v){
   for(u32 i=0; i<v.size(); i++) {
      printf("%d: ",i); v[i].Print(); printf(" \n");
      }
   };
   bool HSameQ(Group &other,u32 start);
   bool CongruenceQ();
	~Group();
};

u32 gpStandardize(u32 mu, u32*S, u32*O, u32*map, u32*Sn, u32*On);
u32 gpFromGens(vector<Matrix> matl, u32*&S, u32*&O);
void gpCosetsGens(u32 mu, u32*S, u32*O, u32*cs, u32*gen);
u32 gpCongruenceTest(u32 n, u32*S, u32*O, u32*rel);

