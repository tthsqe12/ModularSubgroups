u32 GCD(u32 a, u32 b) {
   for(;;){ if(a==0){return b;} b%=a; if(b==0){return a;} a%=b;}
}
u64 GCD(u64 a, u64 b) {
   for(;;){ if(a==0){return b;} b%=a; if(b==0){return a;} a%=b;}
}
u32 LCM(u32 a, u32 b) {
   u32 t=GCD(a,b); return t ? ((a/t)*b) : 0;
}
u64 LCM(u64 a, u64 b) {
   u64 t=GCD(a,b); return t ? ((a/t)*b) : 0;
}
inline i64 Sign(i64 x) {
   return (x != 0) || (x >> 63);
}
inline i64 FlipSign(i64 x, i64 y) {
   return (x ^ (y >> 63)) - (y >> 63);
}
inline i64 Abs(i64 x) {
   return FlipSign(x,x);
}
