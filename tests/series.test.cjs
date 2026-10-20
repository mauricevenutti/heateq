const {test}=require("node:test");
const assert=require("node:assert/strict");
const S=require("../solver.js");
const F=require("../series.js");
function model(source,kind="dirichlet") {
  return S.solve({source,kind,L:Math.PI,N:200,centered:true});
}
test("Parabola has numerical partial sums that update with the Fourier order",()=>{
  const m=model("1-(x/L)^2");
  assert.equal(F.identify(m),null);
  const one=F.formula(m,1),three=F.formula(m,3);
  assert.equal(one.exact,false);assert.match(one.tex,/\\approx/);
  assert.match(one.tex,/1\.03/);
  assert.notEqual(one.tex,three.tex);assert.match(three.tex,/3\\pi/);
});
test("Finite sine sum becomes exact only when all its modes are included",()=>{
  const m=model("sin(pi*(x+L)/(2*L)) + 0.4*sin(3*pi*(x+L)/(2*L))");
  assert.equal(F.formula(m,1).exact,false);
  const full=F.formula(m,3);
  assert.equal(full.exact,true);assert.ok(!full.tex.includes("\\approx"));
  assert.match(full.tex,/0\.40/);
});
test("Rounded coefficients use two decimal places and do not claim exact equality",()=>{
  const fraction=F.formula(model("sin(pi*(x+L)/(2*L))/3"),1);
  assert.equal(fraction.exact,false);assert.match(fraction.tex,/0\.33/);
  const decimal=F.formula(model("0.3333333333333333*sin(pi*(x+L)/(2*L))"),1);
  assert.equal(decimal.exact,false);assert.match(decimal.tex,/0\.33/);
  const negative=F.formula(model("-0.126*sin(pi*(x+L)/(2*L))"),1);
  assert.equal(negative.exact,false);assert.match(negative.tex,/-0\.13/);
});
test("Neumann and periodic exact sums retain their mean",()=>{
  for(const [kind,source] of [
    ["neumann","1+0.6*cos(pi*(x+L)/(2*L))"],
    ["periodic","0.5+sin(pi*x/L)+0.4*cos(2*pi*x/L)"]
  ]) {
    const m=model(source,kind);
    assert.equal(F.formula(m,2).exact,true);
    assert.match(F.formula(m,2).tex,kind==="neumann"?/1\.00 \+/:/0\.50 \+/);
  }
});
test("A small omitted high-frequency term never passes a sampled-equality shortcut",()=>{
  const m=model("sin(pi*(x+L)/(2*L))+1e-14*sin(50*pi*(x+L)/(2*L))");
  assert.equal(F.formula(m,10).exact,false);
  assert.equal(F.formula(m,50).exact,false);
  assert.ok(!F.formula(m,50).tex.includes("0.00"));
});
test("Unshifted centered sine maps to the second eigenmode with the correct sign",()=>{
  const m=model("sin(pi*x/L)"),terms=F.identify(m);
  assert.equal(terms[0].n,2);assert.equal(terms[0].sign,-1);
  assert.equal(F.formula(m,1).exact,false);
  assert.equal(F.formula(m,2).exact,true);assert.match(F.formula(m,2).tex,/=\s*-\\sin/);
});
test("Zero and constant inputs display exact formulas",()=>{
  assert.equal(F.formula(model("0"),1).exact,true);
  assert.equal(F.formula(model("2","neumann"),1).exact,true);
});
