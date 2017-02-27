
void pmId(u32 n, u32*p) {
   u32 i=0; do {p[i]=i;} while (++i<n);
}

void pmMov(u32 n, u32*p, u32*q) {
   u32 i=0; do {p[i]=q[i];} while (++i<n);
}

void pmInv(u32 n, u32*p, u32*q) {
   u32 i=0; do {p[q[i]]=i;} while (++i<n);
}

void pmMulInv(u32 n, u32*p, u32*q, u32*r) {
   u32 i=0; do {p[r[i]]=q[i];} while (++i<n);
}
bool pmMulInvIdQ(u32 n, u32*p, u32*q, u32*r) {
   bool x = true;
   u32 a, b;
   u32 i=0; do {a=r[i]; b=q[i]; p[a]=b; x&=(a==b);} while (++i<n);
   return x;
}


void pmMul(u32 n, u32*p, u32*q, u32*r) {
   u32 i=0; do {p[i]=q[r[i]];} while (++i<n);
}
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s) {
   u32 i=0; do {p[i]=q[r[s[i]]];} while (++i<n);
}
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s, u32*t) {
   u32 i=0; do {p[i]=q[r[s[t[i]]]];} while (++i<n);
}
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s, u32*t, u32*u) {
   u32 i=0; do {p[i]=q[r[s[t[u[i]]]]];} while (++i<n);
}

void pmPow (u32 n, u32*x, u32*y, i32 p, u32*b, u32*t) {
   if (p<0) {p=-p;pmInv(n,b,y);}
   else {pmMov(n,b,y);}
   pmId(n,x);
   while (p>0) {
      if (p&1) {pmMul(n,t,x,b);pmMov(n,x,t);}
      p=p/2;
      pmMul(n,t,b,b);pmMov(n,b,t);
   }
}

bool pmIdQ(u32 n, u32 *p){
   u32 i=0; do {if (p[i]!=i) return false;} while (++i<n);
   return true;
}

bool pmSameQ(u32 n, u32 *p, u32 *q){
   u32 i=0; do {if (p[i]!=q[i]) return false;} while (++i<n);
   return true;
}

u32 pmFixedCount(u32 n, u32*p) {
   u32 ct=0;
   u32 i=0; do {ct+=(p[i]==i);} while (++i<n);
   return ct;
}

u32 pmCycleCount(u32 n, u32*p, u32*vis) {
   u32 i, j, k, ct = 0;
   for (i=0;i<n;i++) {vis[i]=0;};
   for (k=0;k<n;k++) {
      if (vis[k]==0) {
         vis[k]=1;
         for(j=p[k]; j!=k; j=p[j]) {
            assert(vis[j]==0);
            vis[j]=1;
         }
         ct++;
      }
   }
   return ct;
}

u32 pmOrder(u32 n, u32*p, u32*vis) {
   u32 i, j, k, ct = 0, lcm=1;
   for(i=0;i<n;i++){vis[i]=0;};
   for(k=0;k<n;k++) {
      if (vis[k]==0) {
         vis[k]=1;
         ct=1;
         for(j=p[k]; j!=k; j=p[j]) {
            assert(vis[j]==0);
            vis[j]=1;
            ct++;
         }
         lcm = LCM(lcm,ct);
      }
   }
   return lcm;
}

bool pmValidQ(u32 n, u32*p) {
   u32 i, j, k;
   bool * vis = new bool[n];
   for (i=0;i<n;i++){vis[i]=false;};
   for (i=0;i<n;i++){
      if ((p[i]>=n)||vis[p[i]]){ delete[] vis; return false;}
      vis[p[i]]=true;
   };
   delete[] vis;
   return true;
}

void pmPrintCycles(u32 n, u32*p) {
   u32 i, j, k, ct=0;
   bool * visited = new bool[n];
   assert(pmValidQ(n,p));
   //printf("{");
   for(i=0;i<n;i++){visited[i]=false;};
   for(k=0;k<n;k++) {
      if (!visited[k]) {
         printf("(%d",k+1); ct++; visited[k]=true;
         for(j=p[k]; j!=k; j=p[j]) {
            printf(",%d",j+1); ct++; visited[j]=true;
         }
         printf(")");
         //if (ct!=n) {printf(", ");}
      }
   }
   //printf("}");
   delete[] visited;
   return;
}
