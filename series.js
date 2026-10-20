/* Display partial sums. Exact equality is only used for recognized finite
   linear combinations of the actual eigenfunctions, never from a sampled fit. */
(function(root,factory) {
  if(typeof module==="object"&&module.exports)module.exports=factory(require("./vendor/math.js"));
  else root.FourierDisplay=factory(root.math);
})(typeof globalThis!=="undefined"?globalThis:this,function(math) {
  "use strict";
  const unwrap=node=>node.isParenthesisNode?unwrap(node.content):node;
  const canonical=node=>unwrap(node).toString({parenthesis:"auto",implicit:"show"});
  function scalar(node) {
    let result=true;
    node.traverse(part=>{if(part.isSymbolNode&&part.name==="x")result=false;});
    return result;
  }
  function wave(model,n) {
    const multiple=n===1?"":String(n);
    if(model.kind==="periodic")return `\\frac{${multiple}\\pi x}{L}`;
    return model.left===0?`\\frac{${multiple}\\pi x}{L}`:`\\frac{${multiple}\\pi(x+L)}{2L}`;
  }
  function identify(model) {
    const angles=new Map();
    for(let n=1;n<=model.N;n++) {
      const prefix=n===1?"":`${n}*`;
      const source=model.kind==="periodic"||model.left===0
        ?`${prefix}pi*x/L`:`${prefix}pi*(x+L)/(2*L)`;
      angles.set(canonical(math.parse(source)),{n,sign:1});
      // On a centered Dirichlet/Neumann rod, unshifted integer waves are the
      // even modes with phase (-1)^n. This identity is exact, not fitted.
      if(model.kind!=="periodic"&&model.left===-model.L&&2*n<=model.N)
        angles.set(canonical(math.parse(`${prefix}pi*x/L`)),{n:2*n,sign:n%2?-1:1});
    }
    const terms=[];
    function visit(raw,sign=1,factors=[],divisors=[]) {
      const node=unwrap(raw);
      if(scalar(node)) {
        if(node.isConstantNode&&node.value===0)return;
        if(model.kind==="dirichlet")throw new Error("A nonzero constant is not a sine mode.");
        terms.push({n:0,kind:"constant",sign,factors:[...factors,node],divisors});return;
      }
      if(node.isOperatorNode) {
        const [a,b]=node.args;
        if(node.op==="+") {visit(a,sign,factors,divisors);if(b)visit(b,sign,factors,divisors);return;}
        if(node.op==="-") {
          if(!b)visit(a,-sign,factors,divisors);
          else {visit(a,sign,factors,divisors);visit(b,-sign,factors,divisors);}
          return;
        }
        if(node.op==="*"&&(scalar(a)||scalar(b))) {
          const coefficient=scalar(a)?a:b,rest=scalar(a)?b:a;
          if(unwrap(coefficient).isConstantNode&&unwrap(coefficient).value===0)return;
          visit(rest,sign,[...factors,coefficient],divisors);return;
        }
        if(node.op==="/"&&scalar(b)) {visit(a,sign,factors,[...divisors,b]);return;}
      }
      if(node.isFunctionNode&&node.args.length===1) {
        const kind=node.fn.name,mode=angles.get(canonical(node.args[0]));
        if(mode&&((kind==="sin"&&model.kind!=="neumann")||(kind==="cos"&&model.kind!=="dirichlet"))) {
          terms.push({n:mode.n,kind,sign:sign*mode.sign,factors,divisors});return;
        }
      }
      throw new Error("Not a recognized finite Fourier polynomial.");
    }
    try {visit(math.parse(model.source));return terms.sort((a,b)=>a.n-b.n);}
    catch {return null;}
  }
  function exactTerm(model,term) {
    const scope=new Map([["L",model.L]]);
    let value=term.sign;
    for(const node of term.factors)value*=node.compile().evaluate(scope);
    for(const node of term.divisors)value/=node.compile().evaluate(scope);
    const coefficient=numericTex(Math.abs(value));
    const body=term.n===0?coefficient:
      `${coefficient==="1.00"?"":coefficient+"\\,"}\\${term.kind}\\!\\left(${wave(model,term.n)}\\right)`;
    return {sign:value<0?-1:1,body,zero:coefficient==="0.00",
      rounded:Number(value.toFixed(2))!==value};
  }
  function numericTex(value) {
    return value.toFixed(2).replace(/e([+-]?\d+)/,(_,n)=>`\\cdot10^{${Number(n)}}`);
  }
  function formula(model,order,exactTerms=identify(model)) {
    let exact=exactTerms!==null&&exactTerms.every(term=>term.n<=order);
    const terms=[];
    if(exactTerms!==null) {
      const displayed=exactTerms.filter(term=>term.n<=order).map(term=>exactTerm(model,term));
      if(displayed.some(term=>term.rounded))exact=false;
      terms.push(...displayed.filter(term=>!term.zero));
    } else {
      const threshold=model.scale*1e-12;
      function add(value,n,kind) {
        if(Math.abs(value)<=threshold)return;
        const coefficient=numericTex(Math.abs(value));
        if(coefficient==="0.00")return;
        terms.push({sign:value<0?-1:1,body:n===0?coefficient:
          `${coefficient==="1.00"?"":coefficient+"\\,"}\\${kind}\\!\\left(${wave(model,n)}\\right)`});
      }
      if(model.kind!=="dirichlet")add(model.a[0],0,"constant");
      for(let n=1;n<=order;n++) {
        if(model.kind!=="dirichlet")add(model.a[n],n,"cos");
        if(model.kind!=="neumann")add(model.b[n],n,"sin");
      }
    }
    if(!terms.length)terms.push({sign:1,body:"0"});
    const rows=[];
    for(let i=0;i<terms.length;i+=3) {
      const line=terms.slice(i,i+3).map((term,j)=>
        `${term.sign<0?"-":i+j>0?"+":""}${term.body}`).join(" ");
      rows.push((i===0?`\\varphi(x)&${exact?"=":"\\approx"} `:"&\\quad ")+line);
    }
    return {exact,tex:`\\begin{aligned}${rows.join("\\\\[5pt]")}\\end{aligned}`};
  }
  return {identify,formula};
});
