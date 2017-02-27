bool pmIdQ(u32 n, u32 *p);
bool pmSameQ(u32 n, u32 *p, u32 *q);
bool pmValidQ(u32 n, u32*p);
void pmId(u32 n, u32*p);
void pmMov(u32 n, u32*p, u32*q);
void pmPow (u32 n,u32*x, u32*y, i32 p, i32*b, u32*t);
void pmInv(u32 n, u32*p, u32*q);
void pmMulInv(u32 n, u32*p, u32*q, u32*r);
bool pmMulInvIdQ(u32 n, u32*p, u32*q, u32*r);
void pmMul(u32 n, u32*p, u32*q, u32*r);
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s);
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s, u32*t);
void pmMul(u32 n, u32*p, u32*q, u32*r, u32*s, u32*t, u32*u);
void pmPrintCycles(u32 n, u32*p);
u32 pmFixedCount(u32 n, u32*p);
u32 pmCycleCount(u32 n, u32*p, u32*vis);
u32 pmOrder(u32 n, u32*p, u32*vis);


