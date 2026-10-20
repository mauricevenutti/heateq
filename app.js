/* Open index.html directly. All problems use [-L,L], with fixed L = pi. */
(() => {
  "use strict";
  const $ = id => document.getElementById(id);
  const S = window.HeatSolver;
  const LENGTH=Math.PI, REFERENCE_ORDER=200, MAX_ORDER=10, PLOT_SAMPLES=2400;
  const PLAYBACK_MS=30000, ROD_SEGMENTS=240;
  const presets = {
    dirichlet: [
      ["Parabel", "1-(x/L)^2"],
      ["Dreieck", "1-abs(x/L)", "0"],
      ["Wärmepuls", "abs(x) < L/4 ? 1 : 0", "-L/4, L/4"],
      ["Erste Eigenfunktion", "sin(pi*(x+L)/(2*L))"],
      ["Drei Sinusmoden", "sin(pi*(x+L)/(2*L))+0.4*sin(3*pi*(x+L)/(2*L))+0.2*sin(7*pi*(x+L)/(2*L))"]
    ],
    neumann: [
      ["Kosinus + Mittelwert", "1+0.6*cos(pi*(x+L)/(2*L))"],
      ["Glatter Hügel", "(1-(x/L)^2)^2"],
      ["Temperatursprung", "x < 0 ? 1 : 0", "0"],
      ["Wärmepuls", "abs(x) < L/4 ? 1 : 0", "-L/4, L/4"],
      ["Drei Kosinusmoden", "1+0.6*cos(pi*(x+L)/(2*L))+0.3*cos(3*pi*(x+L)/(2*L))+0.2*cos(5*pi*(x+L)/(2*L))"]
    ],
    periodic: [
      ["Sinus", "sin(pi*x/L)"],
      ["Kosinus + Mittelwert", "1+0.6*cos(pi*x/L)"],
      ["Glatter Hügel", "(1-(x/L)^2)^2"],
      ["Wärmepuls", "abs(x) < L/4 ? 1 : 0", "-L/4, L/4"],
      ["Gemischte Moden", "0.5+sin(pi*x/L)+0.4*cos(2*pi*x/L)+0.2*sin(5*pi*x/L)"]
    ]
  };
  let model=null, valid=false, order=0, time=0, maxTime=5, animation=false, animationReady=false;
  let frame=0, previousStamp=0, inputTimer=0, presetBreaks="", initial=[], range=[-1,1];
  let modeSamples=[], referenceTime=-1, referenceCurve=[];
  let temperatureRange=[0,1], rodStops=[], rodEnd=null;
  let exactSeries=null;
  const number=(v,d=5)=>Number(v.toPrecision(d)).toString();
  function tex(target,source,display=false) {
    const el=typeof target==="string"?$(target):target;
    katex.render(source,el,{displayMode:display,throwOnError:true,trust:false,strict:"ignore"});
  }
  function stop() {
    animation=false;cancelAnimationFrame(frame);previousStamp=0;
    $("play").textContent="▶ Start";$("play").setAttribute("aria-label","Animation starten");
  }
  function resetAnimation() {
    stop();time=0;maxTime=5;animationReady=false;
    $("time-slider").value=0;$("time-slider").setAttribute("aria-valuetext","t = 0");
    $("time-value").value=0;$("time-value").setCustomValidity("");
    $("animation-step").hidden=true;
  }
  function resetFourier() {
    resetAnimation();order=0;
    $("terms").value=1;
    for(const id of ["fourier-step","sum-legend","next-animation","reference-note"]) $(id).hidden=true;
    tex("reference-label","\\varphi(x)");
  }
  function invalidateFunction() {
    clearTimeout(inputTimer);resetFourier();valid=false;model=null;
    $("plot-section").hidden=true;$("function-preview").hidden=true;
    $("error").hidden=true;$("expression").removeAttribute("aria-invalid");
  }
  function problem() {
    const kind=$("boundary").value;
    const rb=kind==="dirichlet"?"u(-L,t)=u(L,t)=0":kind==="neumann"?"u_x(-L,t)=u_x(L,t)=0":"u(-L,t)=u(L,t),\\quad u_x(-L,t)=u_x(L,t)";
    tex("problem-math",`\\begin{cases}u_t=u_{xx},&x\\in (-L,L),\\ t>0,\\\\[4pt]${rb},&t>0,\\\\[4pt]u(x,0)=\\varphi(x),&x\\in [-L,L].\\end{cases}`,true);
  }
  function showError(text,formula) {
    $("error").textContent=text;
    if(formula) { const span=document.createElement("span");tex(span,formula);$("error").append(" ",span); }
    $("error").hidden=false;$("expression").setAttribute("aria-invalid","true");
  }
  function sample(fn) {
    return Array.from({length:PLOT_SAMPLES+1},(_,i)=>{const x=-LENGTH+2*LENGTH*i/PLOT_SAMPLES;return {x,y:fn(x)};});
  }
  function prepareCurves() {
    initial=sample(model.f);
    modeSamples=Array.from({length:REFERENCE_ORDER+1},(_,n)=>
      Float64Array.from(initial,p=>S.component(model,n,p.x,0)));
    referenceTime=-1;referenceCurve=[];
    // Choose the viewport once, including every selectable partial sum and the
    // reference at t=0. Heat evolution stays inside these extrema. Never rescale
    // on slider changes or during animation.
    const values=initial.map(p=>p.y);
    const sum=modeSamples[0].slice();
    for(let n=1;n<=REFERENCE_ORDER;n++) {
      for(let i=0;i<sum.length;i++)sum[i]+=modeSamples[n][i];
      if(n<=MAX_ORDER||n===REFERENCE_ORDER) {
        values.push(Math.min(...sum),Math.max(...sum));
      }
    }
    const lo=Math.min(0,...values),hi=Math.max(0,...values),pad=(hi-lo)*.1;
    range=lo===hi?[-1,1]:[lo-pad,hi+pad];
    const temperatures=initial.map(p=>p.y);
    const cold=Math.min(0,...temperatures),hot=Math.max(0,...temperatures);
    temperatureRange=cold===hot?[0,1]:[cold,hot];
    prepareRod();
  }
  function configureTimeline() {
    // The median nonconstant spectral energy selects a relevant decay time.
    // A pure high-frequency input must not inherit the much longer first-mode time.
    const amplitudes=Array.from({length:REFERENCE_ORDER},(_,i)=>
      Math.hypot(model.a[i+1],model.b[i+1]));
    const largest=Math.max(...amplitudes);
    let dominant=1;
    if(largest>Math.max(model.scale,1e-100)*1e-10) {
      const weights=amplitudes.map(a=>(a/largest)**2);
      const halfway=weights.reduce((sum,value)=>sum+value,0)/2;
      let cumulative=0;
      for(let i=0;i<weights.length;i++) {
        cumulative+=weights[i];
        if(cumulative>=halfway){dominant=i+1;break;}
      }
    }
    maxTime=1.5*model.tau/(dominant*dominant);
    $("time-slider").dataset.timeMax=maxTime;
  }
  function temperatureColor(value) {
    const [cold,hot]=temperatureRange;
    const fraction=Math.max(0,Math.min(1,(value-cold)/(hot-cold)));
    return `rgb(${Math.round(224*fraction)},${Math.round(144-123*fraction)},${Math.round(208-175*fraction)})`;
  }
  function prepareRod() {
    $("rod-caption").textContent=model.kind==="periodic"?"Temperatur im Ringstab (abgerollt)":"Temperatur im Stab";
    $("temperature-min").textContent=number(temperatureRange[0],4);
    $("temperature-max").textContent=number(temperatureRange[1],4);
    const stops=Array.from({length:ROD_SEGMENTS+1},(_,i)=>
      `<stop offset="${100*i/ROD_SEGMENTS}%" stop-color="#0090d0"/>`).join("");
    $("rod-plot").innerHTML=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 920 135" role="img" aria-label="Temperatur im Stab" data-temperature-min="${temperatureRange[0]}" data-temperature-max="${temperatureRange[1]}">
      <title>Temperatur der Referenzlösung entlang des Stabs</title>
      <defs>
        <linearGradient id="rod-temperature" gradientUnits="userSpaceOnUse" x1="76" x2="895" y1="0" y2="0">${stops}</linearGradient>
        <linearGradient id="rod-shading" gradientUnits="userSpaceOnUse" x1="0" x2="0" y1="28" y2="92">
          <stop offset="0" stop-color="white" stop-opacity=".28"/>
          <stop offset=".4" stop-color="white" stop-opacity="0"/>
          <stop offset="1" stop-color="#15262e" stop-opacity=".18"/>
        </linearGradient>
        <path id="rod-body" d="M76 28H895A16 32 0 0 1 895 92H76A16 32 0 0 1 76 28Z"/>
      </defs>
      <use href="#rod-body" fill="url(#rod-temperature)"/>
      <use href="#rod-body" fill="url(#rod-shading)"/>
      <use href="#rod-body" fill="none" stroke="#263b41" stroke-width="2.4"/>
      <path d="M76 28A16 32 0 0 1 76 92" fill="none" stroke="#263b41" stroke-width="2" stroke-dasharray="7 6"/>
      <ellipse id="rod-end" cx="895" cy="60" rx="16" ry="32" fill="#0090d0" stroke="#263b41" stroke-width="2.4"/>
      <ellipse cx="895" cy="60" rx="16" ry="32" fill="url(#rod-shading)" pointer-events="none"/>
      <g font-family="Georgia,serif" font-size="14" fill="#36433f" text-anchor="middle">
        <text x="76" y="121">−L</text><text x="485.5" y="121">0</text><text x="895" y="121">L</text>
      </g>
    </svg>`;
    rodStops=Array.from($("rod-plot").querySelectorAll("#rod-temperature stop"));
    rodEnd=$("rod-end");
  }
  function drawRod(curve) {
    for(let i=0;i<rodStops.length;i++) {
      const point=curve[Math.round(i*PLOT_SAMPLES/ROD_SEGMENTS)];
      rodStops[i].setAttribute("stop-color",temperatureColor(point.y));
    }
    rodEnd.setAttribute("fill",temperatureColor(curve[curve.length-1].y));
    const svg=$("rod-plot").firstElementChild;
    svg.setAttribute("aria-label",`Temperatur im Stab bei t = ${number(time,5)}`);
    svg.dataset.time=time;
  }
  function evolvedCurve(n) {
    const values=modeSamples[0].slice();
    for(let k=1;k<=n;k++) {
      const damping=Math.exp(-k*k*time/model.tau);
      for(let i=0;i<values.length;i++)values[i]+=modeSamples[k][i]*damping;
    }
    return initial.map((p,i)=>({x:p.x,y:values[i]}));
  }
  function reference() {
    // The actual input is shown exactly at t=0. For t>0 use the heat series
    // with REFERENCE_ORDER frequencies, independently of the user's slider.
    if(time===0)return initial;
    if(referenceTime!==time) {
      referenceCurve=evolvedCurve(REFERENCE_ORDER);referenceTime=time;
    }
    return referenceCurve;
  }
  function acceptFunction() {
    invalidateFunction();
    if(!$("boundary").value||!$("expression").value.trim())return;
    try {
      model=S.solve({source:$("expression").value,kind:$("boundary").value,L:LENGTH,alpha:1,N:REFERENCE_ORDER,breaks:presetBreaks,centered:true});
      exactSeries=FourierDisplay.identify(model);
      tex("function-preview",`\\varphi(x)=${model.tex}`,true);$("function-preview").hidden=false;
      prepareCurves();$("plot-section").hidden=false;draw();
      if(S.compatibility(model)) {
        const requirement={dirichlet:"\\varphi(-L)=\\varphi(L)=0",neumann:"\\varphi'(-L)=\\varphi'(L)=0",periodic:"\\varphi(-L)=\\varphi(L),\\quad\\varphi'(-L)=\\varphi'(L)"}[model.kind];
        showError("Randbedingungen nicht erfüllt. Erforderlich:",requirement);return;
      }
      valid=true;$("fourier-step").hidden=false;
      configureTimeline();
      chooseOrder();
    } catch(error) { showError(error.message); }
  }
  function chooseOrder() {
    stop();
    order=Number($("terms").value);
    if(!valid||!Number.isInteger(order)||order<1||order>MAX_ORDER)return;
    $("sum-legend").hidden=false;$("next-animation").hidden=animationReady;
    tex("sum-label",`u_{${order}}(x,${animationReady?"t":"0"})`);
    tex("order-value",`n=${order}`);$("terms").setAttribute("aria-valuetext",`${order} Fourierkomponenten`);
    const series=FourierDisplay.formula(model,order,exactSeries);
    tex("series-preview",series.tex,true);
    $("series-preview").dataset.exact=series.exact;
    draw();
  }
  const esc=s=>String(s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&apos;"}[c]));
  function ticks(lo,hi,count) {
    const raw=(hi-lo)/count,magnitude=10**Math.floor(Math.log10(raw));
    const fraction=raw/magnitude,step=(fraction<=1?1:fraction<=2?2:fraction<=2.5?2.5:fraction<=5?5:10)*magnitude;
    const result=[];
    for(let i=Math.ceil(lo/step);i<=Math.floor(hi/step)&&result.length<20;i++)result.push(i*step);
    return result;
  }
  function draw() {
    if(!model)return;
    const width=920,height=400,left=76,right=25,top=25,bottom=55,pw=width-left-right,ph=height-top-bottom;
    const [lo,hi]=range,sx=x=>left+(x+LENGTH)/(2*LENGTH)*pw,sy=y=>top+(hi-y)/(hi-lo)*ph;
    function path(points,jumps=false) {
      let result="",previous=null;
      for(const {x,y} of points) {
        const start=previous===null||(jumps&&Math.abs(y-previous)>(hi-lo)*.3);
        result+=`${start?"M":"L"}${sx(x).toFixed(2)},${sy(y).toFixed(2)}`;previous=y;
      }
      return result;
    }
    let html=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" data-y-min="${lo}" data-y-max="${hi}" role="img" aria-label="${order?`Referenz und Fourierapproximation mit n = ${order}, t = ${number(time)}`:"Anfangsfunktion"}"><title>Temperaturverlauf auf [−L,L]</title><rect width="${width}" height="${height}" fill="white"/><defs><clipPath id="plot-clip"><rect x="${left}" y="${top}" width="${pw}" height="${ph}"/></clipPath></defs><g data-axes font-family="Georgia,serif" font-size="14" fill="#36433f">`;
    for(const y of ticks(lo,hi,5)) {
      const py=sy(y);
      html+=`<path d="M${left} ${py}H${width-right}" stroke="${y===0?"#a9b5ad":"#e5e9e4"}" stroke-width=".8"/><path d="M${left-5} ${py}h5" stroke="#52615a"/><text x="${left-10}" y="${py+5}" text-anchor="end">${esc(number(y,4))}</text>`;
    }
    for(const [x,label] of [[-LENGTH,"−L"],[-LENGTH/2,"−L/2"],[0,"0"],[LENGTH/2,"L/2"],[LENGTH,"L"]]) {
      const px=sx(x);
      html+=`<path d="M${px} ${top}V${top+ph}" stroke="#e5e9e4" stroke-width=".8"/><path d="M${px} ${top+ph}v5" stroke="#52615a"/><text x="${px}" y="${top+ph+23}" text-anchor="middle">${label}</text>`;
    }
    html+=`<path d="M${left} ${top}V${top+ph}H${width-right}" fill="none" stroke="#52615a"/><text x="${left+pw/2}" y="${height-5}" text-anchor="middle" font-style="italic">Position x</text><text transform="translate(20 ${top+ph/2}) rotate(-90)" text-anchor="middle" font-style="italic">Temperatur u</text>`;
    html+='</g>';
    if(order)html+=`<text x="${width-right}" y="16" text-anchor="end" font-family="Georgia,serif" font-size="14">n = ${order} · t = ${number(time,4)}</text>`;
    html+='<g clip-path="url(#plot-clip)">';
    const referenceValues=reference();
    html+=`<path data-series="reference" d="${path(referenceValues,time===0)}" fill="none" stroke="#606a64" stroke-width="2.3" stroke-dasharray="7 5"/>`;
    if(order)html+=`<path data-series="solution" d="${path(evolvedCurve(order))}" fill="none" stroke="#008a94" stroke-width="2.3"/>`;
    html+='</g></svg>';
    $("main-plot").innerHTML=html;
    drawRod(referenceValues);
  }
  function setTime(value,updateInput=true) {
    if(!valid||!order||!Number.isFinite(value)||value<0)return;
    time=value;
    if(time>maxTime){maxTime=time;$("time-slider").dataset.timeMax=maxTime;}
    $("time-slider").value=Math.sqrt(time/maxTime);
    $("time-slider").setAttribute("aria-valuetext",`t = ${number(time,5)}`);
    if(updateInput)$("time-value").value=number(time,7);
    draw();
  }
  $("boundary").addEventListener("change",()=>{
    invalidateFunction();presetBreaks="";$("expression").value="";
    const kind=$("boundary").value;
    $("problem-math").hidden=!kind;$("function-step").hidden=!kind;
    $("preset").replaceChildren(new Option("Eigene Funktion …",""));
    if(!kind)return;
    problem();presets[kind].forEach(([title],i)=>$("preset").append(new Option(title,String(i))));
  });
  $("preset").addEventListener("change",()=>{
    if($("preset").value==="") {$("expression").value="";presetBreaks="";invalidateFunction();$("expression").focus();return;}
    const [,source,breaks=""]=presets[$("boundary").value][Number($("preset").value)];
    $("expression").value=source;presetBreaks=breaks;acceptFunction();
  });
  $("expression").addEventListener("input",()=>{
    presetBreaks="";$("preset").value="";invalidateFunction();
    inputTimer=setTimeout(acceptFunction,300);
  });
  $("function-form").addEventListener("submit",e=>{e.preventDefault();acceptFunction();});
  $("terms").addEventListener("input",chooseOrder);
  $("next-animation").addEventListener("click",()=>{
    if(!valid||!order)return;
    animationReady=true;$("animation-step").hidden=false;$("next-animation").hidden=true;
    tex("reference-label","u_{\\mathrm{ref}}(x,t)");$("reference-note").hidden=false;
    tex("sum-label",`u_{${order}}(x,t)`);$("play").focus();
  });
  $("time-slider").addEventListener("input",()=>{stop();$("time-value").setCustomValidity("");setTime(maxTime*Number($("time-slider").value)**2);});
  $("time-value").addEventListener("input",()=>{
    stop();const value=Number($("time-value").value);
    if(!$("time-value").value||!Number.isFinite(value)||value<0) {
      $("time-value").setCustomValidity("Bitte eine endliche Zeit t ≥ 0 eingeben.");$("time-value").reportValidity();return;
    }
    $("time-value").setCustomValidity("");setTime(value,false);
  });
  $("reset-time").addEventListener("click",()=>{stop();$("time-value").setCustomValidity("");setTime(0);});
  $("play").addEventListener("click",()=>{
    if(!valid||!order||$("animation-step").hidden)return;
    if(animation){stop();return;}
    if(!$("time-value").checkValidity()){$("time-value").reportValidity();return;}
    if(time>=maxTime)setTime(0);
    animation=true;previousStamp=0;$("play").textContent="Ⅱ Pause";$("play").setAttribute("aria-label","Animation pausieren");
    function tick(stamp) {
      if(!animation)return;
      if(previousStamp) {
        const progress=Math.min(1,Math.sqrt(time/maxTime)+Math.min(stamp-previousStamp,100)/PLAYBACK_MS);
        setTime(maxTime*progress*progress);
      }
      previousStamp=stamp;
      if(time>=maxTime)stop();else frame=requestAnimationFrame(tick);
    }
    frame=requestAnimationFrame(tick);
  });
  document.addEventListener("visibilitychange",()=>{if(document.hidden)stop();});
  document.querySelectorAll("[data-tex]").forEach(el=>tex(el,el.dataset.tex));
})();
