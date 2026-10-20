/* Analytical heat semigroup, numerical Fourier coefficients. Browser + Node. */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./vendor/math.js"));
  else root.HeatSolver = factory(root.math);
})(typeof globalThis !== "undefined" ? globalThis : this, function (math) {
  "use strict";
  const functions = new Set(["sin","cos","tan","exp","log","sqrt","abs","sinh","cosh","tanh","asin","acos","atan","floor","ceil","sign","min","max"]);
  const operators = new Set(["+","-","*","/","^","<","<=",">",">=","==","!=","and","or","not"]);
  const types = new Set(["ConstantNode","SymbolNode","OperatorNode","FunctionNode","ParenthesisNode","ConditionalNode"]);

  function expression(source, L, allowX = true) {
    if (!source.trim()) throw new Error("Bitte eine Anfangsfunktion eingeben.");
    if (source.length > 1000) throw new Error("Der Funktionsterm ist zu lang (maximal 1000 Zeichen).");
    let node;
    try { node = math.parse(source); } catch { throw new Error("Der Term ist nicht lesbar. Beispiel: sin(pi*x/L) oder x < L/2 ? 1 : 0."); }
    let count = 0;
    node.traverse((part) => {
      if (++count > 250 || !types.has(part.type)) throw new Error("Erlaubt sind nur reelle Funktionsterme, keine Zuweisungen, Listen oder Objektzugriffe.");
      if (part.isSymbolNode && !new Set(["L","pi","e",...(allowX ? ["x"] : []),...functions]).has(part.name)) throw new Error(`Unbekanntes Symbol „${part.name}“. Erlaubte Variablen: ${allowX ? "x und L" : "L"}.`);
      if (part.isFunctionNode && (!part.fn.isSymbolNode || !functions.has(part.fn.name))) throw new Error("Diese Funktion ist nicht zugelassen.");
      if (part.isOperatorNode && !operators.has(part.op)) throw new Error(`Der Operator „${part.op}“ wird nicht unterstützt.`);
      if (part.isConstantNode && typeof part.value !== "number") throw new Error("Bitte nur numerische Konstanten verwenden.");
    });
    const compiled = node.compile();
    const scope = new Map([["L",L],["x",0]]);
    const evaluate = (x) => {
      scope.set("x",x);
      let result;
      try { result = compiled.evaluate(scope); } catch { throw new Error(`Die Funktion ist bei x = ${x.toPrecision(6)} nicht auswertbar.`); }
      if (typeof result !== "number" || !Number.isFinite(result) || Math.abs(result) > 1e100) throw new Error(`Die Funktion muss reell und endlich sein (Problem bei x = ${x.toPrecision(6)}).`);
      return result;
    };
    return { evaluate, tex:node.toTex({parenthesis:"auto", implicit:"hide"}), source };
  }
  function breakpoints(source, L, left) {
    if (!source.trim()) return [];
    const values = source.split(",").map(s => expression(s.trim(),L,false).evaluate(0));
    if (values.length > 32) throw new Error("Bitte höchstens 32 Integrationsstellen angeben.");
    if (values.some(v => v <= left || v >= L)) throw new Error("Alle Integrationsstellen müssen im Inneren des räumlichen Intervalls liegen.");
    return [...new Set(values)].sort((a,b)=>a-b);
  }
  // Composite 8-point Gauss–Legendre integration; open nodes avoid endpoint artefacts.
  const gx = [-.9602898564975363,-.7966664774136267,-.525532409916329,-.1834346424956498,.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363];
  const gw = [.1012285362903763,.2223810344533745,.3137066458778873,.362683783378362,.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763];
  function integrate(f, kind, L, N, points, panels, left) {
    const a = new Float64Array(N+1), b = new Float64Array(N+1);
    const basisLength = kind === "periodic" ? L : L-left;
    const origin = kind === "periodic" ? 0 : left;
    const bounds = [left,...points,L];
    let samples = 0;
    for (let segment=0;segment<bounds.length-1;segment++) {
      const lo=bounds[segment], hi=bounds[segment+1];
      const cells=Math.max(1,Math.ceil(panels*(hi-lo)/(L-left)));
      const h=(hi-lo)/cells;
      for(let j=0;j<cells;j++) for(let q=0;q<8;q++) {
        const x=lo+(j+.5)*h+gx[q]*h/2;
        // Work in dimensionless integration weights so extreme L does not inflate sums.
        const y=f(x)*gw[q]*(h/basisLength)/2;
        samples++;
        a[0]+=y/(kind==="periodic"?2:1);
        const angle=Math.PI*(x-origin)/basisLength, c=Math.cos(angle), s=Math.sin(angle);
        let cn=1,sn=0;
        for(let n=1;n<=N;n++) {
          const next=cn*c-sn*s;
          sn=sn*c+cn*s; cn=next;
          if(kind!=="dirichlet") a[n]+=y*cn*(kind==="periodic"?1:2);
          if(kind!=="neumann") b[n]+=y*sn*(kind==="periodic"?1:2);
        }
      }
    }
    if(kind==="dirichlet") a[0]=0;
    return {a,b,samples};
  }
  function solve({source,kind="dirichlet",L=1,alpha=1,N=30,breaks="",centered=false}) {
    if (!["dirichlet","neumann","periodic"].includes(kind)) throw new Error("Unbekannte Randbedingung.");
    if(!Number.isFinite(L)||L<1e-6||L>1e6) throw new Error("L muss zwischen 0.000001 und 1000000 liegen.");
    if(!Number.isFinite(alpha)||alpha<1e-6||alpha>1e6) throw new Error("α muss zwischen 0.000001 und 1000000 liegen.");
    if(!Number.isInteger(N)||N<1||N>200) throw new Error("Die Fourierordnung N muss eine ganze Zahl zwischen 1 und 200 sein.");
    const parsed=expression(source,L), left=kind==="periodic"||centered?-L:0;
    const basisLength=kind==="periodic"?L:L-left;
    const origin=kind==="periodic"?0:left;
    const points=breakpoints(breaks,L,left);
    // Validate endpoints as well: this app accepts finite, bounded input functions.
    const f=parsed.evaluate;
    let scale=0;
    for(let i=0;i<=1024;i++) scale=Math.max(scale,Math.abs(f(left+(L-left)*i/1024)));
    let panels=Math.max(32,Math.ceil(N/2));
    let previous=integrate(f,kind,L,N,points,panels,left), current, error=Infinity, converged=false;
    let samples=previous.samples;
    for(let pass=0;pass<6;pass++) {
      panels*=2;
      current=integrate(f,kind,L,N,points,panels,left); samples+=current.samples;
      error=0;
      for(let n=0;n<=N;n++) error=Math.max(error,Math.abs(current.a[n]-previous.a[n]),Math.abs(current.b[n]-previous.b[n]));
      const coeffScale=Math.max(...current.a.map(Math.abs),...current.b.map(Math.abs));
      scale=Math.max(scale,coeffScale);
      if(error <= 1e-9*Math.max(scale,1e-12)) { converged=true; break; }
      previous=current;
    }
    const tau=(basisLength/(Math.PI*alpha))**2;
    return {source,kind,L,alpha,N,left,origin,basisLength,points,f,tex:parsed.tex,...current,error,converged,samples,scale,tau};
  }
  function component(model,n,x,t) {
    if(n===0) return model.a[0];
    const angle=n*Math.PI*(x-model.origin)/model.basisLength;
    return Math.exp(-n*n*t/model.tau)*(model.a[n]*Math.cos(angle)+model.b[n]*Math.sin(angle));
  }
  function value(model,x,t,order=model.N) {
    let sum=model.a[0];
    for(let n=1;n<=order;n++) sum+=component(model,n,x,t);
    return sum;
  }
  function compatibility(model) {
    const {f,left,L,kind}=model, span=L-left, h=span*1e-5;
    const tolerance=1e-5*Math.max(model.scale,1e-12);
    const d0=(-3*f(left)+4*f(left+h)-f(left+2*h))/(2*h);
    const dL=(3*f(L)-4*f(L-h)+f(L-2*h))/(2*h);
    if(kind==="dirichlet" && (Math.abs(f(left))>tolerance||Math.abs(f(L))>tolerance))
      return "Die Anfangsfunktion hat nichtverschwindende Randwerte. Die Sinusreihe setzt beide Endtemperaturen auf null. Bei t = 0 stimmt sie dort nicht mit der Eingabe überein; für t > 0 gilt die Randbedingung exakt.";
    if(kind==="neumann" && Math.max(Math.abs(d0),Math.abs(dL))*span>tolerance)
      return "Die numerisch geschätzten Randableitungen der Anfangsfunktion sind nicht null. Die Kosinusreihe erfüllt für t > 0 isolierte Enden; eine bis t = 0 verträgliche Randableitung ist hier nicht gegeben.";
    if(kind==="periodic" && (Math.abs(f(left)-f(L))>tolerance || Math.abs(d0-dL)*span>tolerance))
      return "Werte oder numerisch geschätzte Ableitungen passen an der Naht des Ringstabs nicht zusammen. Die Reihe verwendet die periodische Fortsetzung. Bei einem Wertsprung konvergiert sie dort zum Mittelwert der einseitigen Grenzwerte.";
    return "";
  }
  return {expression,breakpoints,solve,component,value,compatibility};
});
