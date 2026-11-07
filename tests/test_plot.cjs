const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const ids={};let strokes=0;
const ctx=new Proxy({}, {get:(o,k)=>o[k]??(()=>{if(k==='stroke')strokes++}),set:(o,k,v)=>(o[k]=v,true)});
class Element {
 constructor(){this.clientWidth=1000;this.clientHeight=370;this.style={};this.children=[];this.hidden=true;this.offsetWidth=220;this.offsetHeight=110;}
 getContext(){return ctx}
 getBoundingClientRect(){return {left:0,top:0,width:this.clientWidth,height:this.clientHeight}}
 addEventListener(type,handler,options){this['on'+type]=handler;this.listenerOptions=options}
 setPointerCapture(){} append(e){this.children.push(e)} replaceChildren(){this.children=[]}
}
let frame;
const sandbox=vm.createContext({document:{getElementById:id=>ids[id]??=new Element(),createElement:()=>new Element()},window:{devicePixelRatio:2},ResizeObserver:class{observe(){}},requestAnimationFrame:f=>(frame=f,1),Blob,URL,setTimeout});
vm.runInContext(fs.readFileSync(__dirname+'/../viewer/plot.js','utf8')+'\nglobalThis.TrackPlot=TrackPlot',sandbox);
let change;const plot=new sandbox.TrackPlot((...args)=>change=args);
plot.setRange(100,200,1000);
plot.navigate(-100,100);assert.equal(JSON.stringify(change.slice(0,2)),'[0,200]');
plot.navigate(950,1050);assert.equal(JSON.stringify(change.slice(0,2)),'[900,1000]');
plot.setRange(100,200,1000);plot.canvas.onwheel({preventDefault(){},clientX:526,deltaY:-100});
assert.ok(change[0]>100&&change[1]<200);assert.equal(change[2],true);
assert.equal(plot.canvas.listenerOptions.passive,false);
// Zoom limits must not trigger useless data requests.
plot.setRange(0,1000,1000);change=null;
plot.canvas.onwheel({preventDefault(){},clientX:526,deltaY:100});assert.equal(change,null);
// Line-based wheel devices zoom around the pointer, just like pixel-based ones.
plot.setRange(100,900,1000);
plot.canvas.onwheel({preventDefault(){},clientX:72,deltaY:-3,deltaMode:1});
assert.equal(change[0],100);assert.ok(change[1]<900);
plot.canvas.ondblclick();assert.equal(JSON.stringify(change.slice(0,2)),'[0,1000]');
plot.setRange(90000,91000,100000);plot.canvas.ondblclick();assert.equal(JSON.stringify(change.slice(0,2)),'[0,100000]');
plot.setRange(0,4,4);const data={start:0,end:4,bin_size:1,bin_start:[0,1,2,3],bin_end:[1,2,3,4],plus:[0,null,5,0],minus:[2,null,0,1],eligible:[1,0,1,1]};
assert.equal(plot.logScale,true);
plot.setData(data);plot.setOverview(data);
let g=plot.geometry(1000,370);assert.equal(g.y(0),g.bottom);assert.ok(Number.isFinite(g.y(0)));
const logFraction=(g.bottom-g.y(1))/(g.bottom-g.y(5));
plot.logScale=false;g=plot.geometry(1000,370);
const linearFraction=(g.bottom-g.y(1))/(g.bottom-g.y(5));
assert.ok(logFraction>linearFraction);assert.ok(Math.abs(linearFraction-.2)<1e-10);
plot.logScale=true;assert.equal(JSON.stringify(plot.ticks(100)),'[0,1,10,100]');
frame();assert.ok(strokes>0);assert.equal(plot.canvas.width,2000);
plot.hover=.65;plot.showTooltip({clientX:600,clientY:100});assert.equal(ids.tooltip.hidden,false);assert.equal(ids.tooltip.children[0].textContent,'Position 2');
plot.pan=true;plot.setRange(20,40,100);plot.canvas.onpointerdown({button:0,pointerId:1,clientX:526});plot.canvas.onpointermove({clientX:753});plot.canvas.onpointerup({});assert.equal(JSON.stringify(change.slice(0,2)),'[15,35]');
// Letters use reference coordinates, including excluded bases, only at readable zoom.
plot.setData({...data,dna:'ACGT'});plot.setRange(0,4,1000);
let letters=plot.nucleotides(plot.geometry(1000,370));
assert.equal(letters.map(letter=>letter.base).join(''),'ACGT');
assert.equal(letters[2].x,plot.geometry(1000,370).x(2));
assert.equal(plot.geometry(1000,370).top,48);
plot.setRange(1,3,1000);assert.equal(plot.nucleotides(plot.geometry(1000,370)).map(l=>l.base).join(''),'CG');
plot.setRange(0,1000,1000);assert.equal(plot.nucleotides(plot.geometry(1000,370)).length,0);
plot.setRange(10,14,1000);assert.equal(plot.nucleotides(plot.geometry(1000,370)).length,0);
plot.setData({...data,dna:null});plot.setRange(0,4,1000);
assert.equal(plot.nucleotides(plot.geometry(1000,370)).length,0);
// Shift selects the alternate action for the whole gesture, without changing the default.
plot.pan=false;plot.setRange(20,40,100);
plot.canvas.onpointerdown({button:0,pointerId:1,clientX:526,shiftKey:true});
plot.canvas.onpointermove({clientX:753,shiftKey:false});
plot.canvas.onpointerup({shiftKey:false});
assert.equal(JSON.stringify(change.slice(0,2)),'[15,35]');assert.equal(plot.pan,false);
// With pan selected, Shift-drag zooms instead.
plot.pan=true;plot.setRange(0,100,100);
plot.canvas.onpointerdown({button:0,pointerId:1,clientX:298.5,shiftKey:true});
plot.canvas.onpointermove({clientX:751.5,shiftKey:true});
plot.canvas.onpointerup({shiftKey:false});
assert.equal(JSON.stringify(change.slice(0,2)),'[25,75]');assert.equal(plot.pan,true);
console.log('Canvas scroll zoom, zoom limits, log/linear scales, pan, reset, hover, and drawing checks passed');
// Scroll bursts update coordinates immediately but issue only one pair of requests.
let viewportChange,fetches=[],timers=new Map(),timerId=0;
const ui=new Proxy({}, {get:(o,id)=>o[id]??=({value:'',clientWidth:1000,parentElement:{setAttribute(){}},setAttribute(){}})});
const viewer=vm.createContext({document:{getElementById:id=>ui[id]},TrackPlot:class{constructor(callback){viewportChange=callback}},AbortController,URLSearchParams,
 fetch:(path,options)=>{fetches.push({path,options});return new Promise(()=>{})},
 setTimeout:callback=>(timers.set(++timerId,callback),timerId),clearTimeout:id=>timers.delete(id)});
vm.runInContext(fs.readFileSync(__dirname+'/../viewer/viewer.js','utf8'),viewer);
vm.runInContext("state.contig={id:'1',name:'chr1'};state.dataset='test';pending.set('view',new AbortController());globalThis.previousView=pending.get('view')",viewer);
viewportChange(100,900,true);viewportChange(200,800,true);
assert.equal(ui.start.value,200);assert.equal(ui.end.value,800);
assert.equal(viewer.previousView.signal.aborted,true);
assert.equal(fetches.length,2);assert.equal(timers.size,1);
[...timers.values()][0]();assert.equal(fetches.length,4);
assert.ok(fetches.slice(2).every(r=>r.path.includes('start=200&end=800')));
console.log('Scroll request coalescing and stale-request cancellation checks passed');
