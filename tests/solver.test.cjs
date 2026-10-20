const {test}=require("node:test");
const assert=require("node:assert/strict");
const S=require("../solver.js");
function near(actual,expected,tol=1e-9) {
  assert.ok(Math.abs(actual-expected)<=tol,`${actual} ≠ ${expected} (tolerance ${tol})`);
}
test("Dirichlet: single sine mode with non-unit L and alpha",()=>{
  const m=S.solve({source:"3*sin(2*pi*x/L)",L:2.7,alpha:.4,N:15});
  near(m.b[2],3);
  for(let n=1;n<=15;n++)if(n!==2)near(m.b[n],0);
  const x=.61,t=.27;
  near(S.value(m,x,t),3*Math.sin(2*Math.PI*x/2.7)*Math.exp(-(.4**2)*(2*Math.PI/2.7)**2*t));
  near(S.value(m,0,t),0);near(S.value(m,2.7,t),0);
  assert.equal(m.converged,true);
});
test("Dirichlet: parabola coefficients 32/(pi^3 n^3) for odd n",()=>{
  const m=S.solve({source:"4*x*(L-x)/L^2",N:30});
  for(let n=1;n<=m.N;n++)near(m.b[n],n%2?32/(Math.PI*n)**3:0);
});
test("Neumann: cosine modes, mean normalization and conservation",()=>{
  const m=S.solve({source:"2 + 0.6*cos(pi*x/L) - 0.2*cos(3*pi*x/L)",kind:"neumann",L:3,alpha:2,N:12});
  near(m.a[0],2);near(m.a[1],.6);near(m.a[3],-.2);
  near(S.value(m,.7,100),2);
  for(let n=0;n<=m.N;n++)near(m.b[n],0);
  for(const t of [0,.03,.4]) {
    let mean=0;
    for(let i=0;i<1000;i++)mean+=S.value(m,3*(i+.5)/1000,t)/1000;
    near(mean,2);
  }
});
test("Periodic: independent sine and cosine coefficients on [-L,L]",()=>{
  const m=S.solve({source:"0.5+sin(pi*x/L)+0.4*cos(2*pi*x/L)+0.2*sin(5*pi*x/L)",kind:"periodic",L:2,alpha:.5,N:12});
  near(m.a[0],.5);near(m.a[2],.4);near(m.b[1],1);near(m.b[5],.2);
  near(S.value(m,-2,.08),S.value(m,2,.08));
  near(S.value(m,.32,0),m.f(.32));
});
test("Discontinuous step: Fourier integrals split at a known discontinuity",()=>{
  const m=S.solve({source:"x < L/3 ? 1 : 0",N:30,breaks:"L/3"});
  for(let n=1;n<=m.N;n++)near(m.b[n],2*(1-Math.cos(n*Math.PI/3))/(n*Math.PI));
  assert.equal(m.converged,true);
});
test("Periodic sawtooth agrees with the script's Fourier coefficients",()=>{
  const m=S.solve({source:"x",kind:"periodic",L:Math.PI,N:30});
  near(m.a[0],0);
  for(let n=1;n<=m.N;n++) {near(m.a[n],0);near(m.b[n],2*(-1)**(n+1)/n);}
  assert.ok(S.compatibility(m).includes("Naht"));
});
test("Finite Fourier solution satisfies the heat PDE",()=>{
  const m=S.solve({source:"x*(L-x)",N:20,L:2,alpha:.7});
  const x=.7,t=.2,h=1e-4;
  const ut=(S.value(m,x,t+h)-S.value(m,x,t-h))/(2*h);
  const uxx=(S.value(m,x+h,t)-2*S.value(m,x,t)+S.value(m,x-h,t))/(h*h);
  near(ut,m.alpha*m.alpha*uxx,1e-6);
});
test("Zero and very small inputs remain well defined",()=>{
  const zero=S.solve({source:"0",N:10});near(S.value(zero,.3,1),0);
  const tiny=S.solve({source:"1e-15*sin(pi*x/L)",N:10});
  near(tiny.b[1],1e-15,1e-25);
});
test("Expression whitelist accepts math and rejects assignments and access",()=>{
  near(S.expression("2x + abs(-3) + sqrt(4)",1).evaluate(.5),6);
  near(S.expression("x > 0 and x < L ? 1 : 0",1).evaluate(.5),1);
  for(const source of ["x=1","import('x')","[1,2]","x.constructor","random()","foo(x)","sin(x); cos(x)","'hello'"]) {
    assert.throws(()=>S.expression(source,1),source);
  }
  for(const source of ["sqrt(-1)","1/(x-.5)","log(0)"])assert.throws(()=>S.solve({source,N:4}),source);
});
test("Rejects invalid model parameters and integration breakpoints",()=>{
  for(const params of [{L:0},{alpha:-1},{N:2.5},{N:201},{L:Infinity},{kind:"unknown"},{breaks:"2*L"},{breaks:"x/2"}]) {
    assert.throws(()=>S.solve({source:"x",...params}));
  }
});
test("Centered Dirichlet rod: translated eigenbasis and correct decay length",()=>{
  const L=Math.PI;
  const m=S.solve({source:"sin(pi*(x+L)/(2*L)) + 0.3*sin(3*pi*(x+L)/(2*L))",L,N:200,centered:true});
  near(m.left,-L);near(m.basisLength,2*L);near(m.tau,4);
  near(m.b[1],1);near(m.b[3],.3);
  for(const x of [-L,-.7,0,1.5,L]) {
    const t=.4;
    const expected=Math.sin((x+L)/2)*Math.exp(-t/4)+.3*Math.sin(3*(x+L)/2)*Math.exp(-9*t/4);
    near(S.value(m,x,t),expected);
  }
  near(S.value(m,-L,.4),0);near(S.value(m,L,.4),0);
  assert.equal(S.compatibility(m),"");
});
test("Centered Neumann rod: mean, coefficients and boundary compatibility",()=>{
  const m=S.solve({source:"2 + 0.6*cos(pi*(x+L)/(2*L))",kind:"neumann",L:Math.PI,N:200,centered:true});
  near(m.a[0],2);near(m.a[1],.6);near(m.tau,4);
  near(S.value(m,.8,.5),2+.6*Math.cos((.8+Math.PI)/2)*Math.exp(-.5/4));
  assert.equal(S.compatibility(m),"");
  const bad=S.solve({source:"x",kind:"neumann",L:Math.PI,N:10,centered:true});
  assert.notEqual(S.compatibility(bad),"");
});
test("Centered periodic rod retains its full-period eigenvalues",()=>{
  const m=S.solve({source:"0.5+sin(pi*x/L)+0.4*cos(2*pi*x/L)",kind:"periodic",L:Math.PI,N:200,centered:true});
  near(m.tau,1);near(m.a[0],.5);near(m.b[1],1);near(m.a[2],.4);
  near(S.value(m,.7,.4),.5+Math.sin(.7)*Math.exp(-.4)+.4*Math.cos(1.4)*Math.exp(-1.6));
  near(S.value(m,-Math.PI,.4),S.value(m,Math.PI,.4));
});
