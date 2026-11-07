/* Bounded-resolution canvas renderer. No per-base DOM or whole-genome arrays. */
class TrackPlot {
  constructor(change) {
    this.canvas=document.getElementById('plot');this.overview=document.getElementById('overview');
    this.tooltip=document.getElementById('tooltip');this.change=change;this.visible={plus:true,minus:true};
    this.range=[0,1];this.length=1;this.logScale=true;this.pan=false;this.drag=null;this.hover=null;
    this.colors={plus:'#2563eb',minus:'#e36a45'};
    new ResizeObserver(()=>this.schedule()).observe(this.canvas);
    new ResizeObserver(()=>this.schedule()).observe(this.overview);
    this.canvas.addEventListener('wheel',event=>{
      event.preventDefault();
      const units=event.deltaMode===1?16:event.deltaMode===2?this.canvas.clientHeight:1;
      const factor=Math.exp(Math.max(-.7,Math.min(.7,event.deltaY*units*.002)));
      const fraction=this.fraction(event,this.canvas),width=this.range[1]-this.range[0];
      const anchor=this.range[0]+width*fraction,newWidth=width*factor;
      this.hover=null;this.tooltip.hidden=true;
      this.navigate(anchor-newWidth*fraction,anchor+newWidth*(1-fraction),true);
    },{passive:false});
    this.canvas.ondblclick=()=>this.navigate(0,this.length);
    this.canvas.onpointerdown=e=>{if(e.button!==0)return;this.canvas.setPointerCapture(e.pointerId);this.drag={x:this.fraction(e,this.canvas),range:[...this.range],pan:this.pan!==!!e.shiftKey};this.canvas.style.cursor=this.drag.pan?'grabbing':'crosshair';this.tooltip.hidden=true};
    this.canvas.onpointermove=e=>{
      if(this.drag){this.drag.to=this.fraction(e,this.canvas);if(this.drag.pan){const width=this.drag.range[1]-this.drag.range[0],delta=(this.drag.x-this.drag.to)*width;this.preview(this.drag.range[0]+delta,this.drag.range[1]+delta)}this.schedule();return}
      this.canvas.style.cursor=(this.pan!==!!e.shiftKey)?'grab':'crosshair';
      this.hover=this.fraction(e,this.canvas);this.schedule();this.showTooltip(e);
    };
    this.canvas.onpointerup=e=>{if(!this.drag)return;const d=this.drag;this.drag=null;this.canvas.style.cursor=(this.pan!==!!e.shiftKey)?'grab':'crosshair';if(d.pan)this.navigate(...this.range);else if(d.to!==undefined&&Math.abs(d.to-d.x)>.004){const width=d.range[1]-d.range[0];this.navigate(d.range[0]+Math.min(d.x,d.to)*width,d.range[0]+Math.max(d.x,d.to)*width)}this.schedule()};
    this.canvas.onpointercancel=()=>{this.drag=null;this.canvas.style.cursor=this.pan?'grab':'crosshair';this.schedule()};
    this.canvas.onpointerleave=()=>{this.hover=null;this.tooltip.hidden=true;this.schedule()};
    this.canvas.onkeydown=e=>{if(['ArrowLeft','ArrowRight','+','=','-'].includes(e.key)){e.preventDefault();if(e.key.startsWith('Arrow'))this.shift(e.key==='ArrowLeft'?-.5:.5);else this.zoom(e.key==='-'?2:.5)}};
    this.overview.onpointerdown=e=>{if(e.button!==0)return;this.overview.setPointerCapture(e.pointerId);const p=this.fraction(e,this.overview,true)*this.length;const tolerance=this.length/this.overview.clientWidth*12;let mode=Math.abs(p-this.range[0])<tolerance?'left':Math.abs(p-this.range[1])<tolerance?'right':'move';this.overDrag={p,range:[...this.range],mode};if(mode==='move'&&(p<this.range[0]||p>this.range[1])){const w=this.range[1]-this.range[0];this.preview(p-w/2,p+w/2);this.overDrag={p,range:[...this.range],mode}}};
    this.overview.onpointermove=e=>{if(!this.overDrag)return;const d=this.overDrag,p=this.fraction(e,this.overview,true)*this.length,delta=p-d.p;if(d.mode==='left')this.preview(Math.min(p,d.range[1]-1),d.range[1]);else if(d.mode==='right')this.preview(d.range[0],Math.max(p,d.range[0]+1));else this.preview(d.range[0]+delta,d.range[1]+delta)};
    this.overview.onpointerup=()=>{if(this.overDrag){this.overDrag=null;this.navigate(...this.range)}};
    this.overview.onpointercancel=()=>{this.overDrag=null;this.schedule()};
  }
  fraction(event,canvas,overview=false){const r=canvas.getBoundingClientRect(),left=overview?0:72,right=overview?0:22;return Math.max(0,Math.min(1,(event.clientX-r.left-left)/(r.width-left-right)))}
  clamp(a,b){const width=Math.min(this.length,Math.max(1,Math.round(b-a)));a=Math.max(0,Math.min(this.length-width,Math.round(a)));return [a,a+width]}
  preview(a,b){this.range=this.clamp(a,b);this.schedule()}
  navigate(a,b,debounce=false){
    const next=this.clamp(a,b);
    if(next[0]===this.range[0]&&next[1]===this.range[1]&&debounce)return;
    this.preview(...next);this.change(...this.range,debounce);
  }
  zoom(factor){const center=(this.range[0]+this.range[1])/2,width=(this.range[1]-this.range[0])*factor;this.navigate(center-width/2,center+width/2)}
  shift(fraction){const delta=(this.range[1]-this.range[0])*fraction;this.navigate(this.range[0]+delta,this.range[1]+delta)}
  setRange(start,end,length){this.length=length;this.range=[start,end];this.hover=null;this.tooltip.hidden=true;this.schedule()}
  setData(data){this.data=data;this.schedule()}
  setOverview(data){this.overData=data;this.schedule()}
  schedule(){if(this.frame)return;this.frame=requestAnimationFrame(()=>{this.frame=null;this.draw();this.drawOverview()})}
  context(canvas){const ratio=window.devicePixelRatio||1,width=canvas.clientWidth,height=canvas.clientHeight;if(canvas.width!==Math.round(width*ratio)||canvas.height!==Math.round(height*ratio)){canvas.width=Math.round(width*ratio);canvas.height=Math.round(height*ratio)}const ctx=canvas.getContext('2d');ctx.setTransform(ratio,0,0,ratio,0,0);ctx.clearRect(0,0,width,height);return [ctx,width,height]}
  scale(value){return this.logScale?Math.log10(1+value):value}
  ticks(peak){
    if(!this.logScale)return Array.from({length:5},(_,i)=>peak*i/4);
    const ticks=[0],step=Math.max(1,Math.ceil(Math.log10(peak+1)/6));
    for(let exponent=0;10**exponent<=peak;exponent+=step)ticks.push(10**exponent);
    return ticks;
  }
  geometry(width,height){const left=72,right=width-22,top=this.data?.dna&&(right-left)/(this.range[1]-this.range[0])>=12?48:22,bottom=height-45;const span=this.range[1]-this.range[0];let peak=1;const d=this.data;if(d)for(let i=0;i<d.plus.length;i++)if(d.bin_end[i]>this.range[0]&&d.bin_start[i]<this.range[1])for(const key of ['plus','minus'])if(this.visible[key])peak=Math.max(peak,d[key][i]||0);return {left,right,top,bottom,peak,x:p=>left+(p-this.range[0])/span*(right-left),y:v=>bottom-this.scale(v)/this.scale(peak)*(bottom-top-14)}}
  nucleotides(g){
    const d=this.data;
    if(!d?.dna||d.bin_size!==1||(g.right-g.left)/(this.range[1]-this.range[0])<12)return [];
    const colors={A:'#18814a',C:'#2563eb',G:'#a96a08',T:'#c43c3c'};
    const letters=[];
    for(let p=Math.max(this.range[0],d.start);p<Math.min(this.range[1],d.start+d.dna.length);p++){
      const base=d.dna[p-d.start];letters.push({base,x:g.x(p),color:colors[base]||'#6b7280'});
    }
    return letters;
  }
  paths(ctx,data,g){
    for(const key of ['plus','minus']){if(!this.visible[key])continue;ctx.beginPath();ctx.strokeStyle=this.colors[key];ctx.lineWidth=1.6;let connected=false;
      for(let i=0;i<data[key].length;i++){
        if(data.bin_end[i]<=this.range[0]||data.bin_start[i]>=this.range[1]){connected=false;continue}
        const value=data[key][i];if(value===null){connected=false;continue}
        const x=g.x((data.bin_start[i]+data.bin_end[i]-1)/2),y=g.y(value);
        if(connected)ctx.lineTo(x,y);else ctx.moveTo(x,y);connected=true;
        if((i===0||data[key][i-1]===null)&&(i===data[key].length-1||data[key][i+1]===null)){ctx.moveTo(x-2,y);ctx.lineTo(x+2,y)}
      }ctx.stroke();
    }
  }
  draw(){const [ctx,w,h]=this.context(this.canvas);if(w<100)return;const g=this.geometry(w,h);this.geo=g;ctx.fillStyle='#ffffff';ctx.fillRect(0,0,w,h);ctx.font='11px system-ui';ctx.textBaseline='middle';
    for(const value of this.ticks(g.peak)){const y=g.y(value);ctx.strokeStyle='#edf1f6';ctx.beginPath();ctx.moveTo(g.left,y);ctx.lineTo(g.right,y);ctx.stroke();ctx.textAlign='right';ctx.fillStyle='#72839a';ctx.fillText(value.toLocaleString(undefined,{maximumFractionDigits:1,notation:value>=1e6?'compact':'standard'}),g.left-12,y)}
    const span=this.range[1]-this.range[0],step=Math.max(1,10**Math.floor(Math.log10(span/5)));const mult=span/5/step;const tick=step*(mult>5?10:mult>2?5:mult>1?2:1);ctx.textAlign='center';
    for(let p=Math.ceil(this.range[0]/tick)*tick;p<this.range[1];p+=tick)ctx.fillText(p.toLocaleString(),g.x(p),g.bottom+20);
    ctx.fillStyle='#60758e';ctx.textAlign='left';ctx.fillText(this.logScale?'csRNA counts · log₁₀(count + 1)':'csRNA counts',g.left,10);ctx.textAlign='center';ctx.fillText('Reference position · 0-based',w/2,h-6);
    const letters=this.nucleotides(g);
    if(letters.length){
      ctx.font='bold 12px ui-monospace, monospace';ctx.textAlign='center';
      for(const letter of letters){ctx.fillStyle=letter.color;ctx.fillText(letter.base,letter.x,31)}
    }
    ctx.save();ctx.beginPath();ctx.rect(g.left,g.top,g.right-g.left,g.bottom-g.top);ctx.clip();
    const d=this.data;if(d){ctx.fillStyle='#f1f3f6';for(let i=0;i<d.eligible.length;i++)if(!d.eligible[i]){const x=g.x(d.bin_start[i]);ctx.fillRect(x,g.top,Math.max(1,g.x(d.bin_end[i])-x),g.bottom-g.top)}this.paths(ctx,d,g)}
    if(this.hover!==null){const x=g.left+this.hover*(g.right-g.left);ctx.strokeStyle='#9baabd';ctx.setLineDash([3,4]);ctx.beginPath();ctx.moveTo(x,g.top);ctx.lineTo(x,g.bottom);ctx.stroke();ctx.setLineDash([])}
    if(this.drag&&!this.drag.pan&&this.drag.to!==undefined){ctx.fillStyle='#2563eb20';ctx.fillRect(g.left+Math.min(this.drag.x,this.drag.to)*(g.right-g.left),g.top,Math.abs(this.drag.to-this.drag.x)*(g.right-g.left),g.bottom-g.top)}ctx.restore();
  }
  drawOverview(){const [ctx,w,h]=this.context(this.overview),d=this.overData;if(!d)return;let peak=1;for(const key of ['plus','minus'])if(this.visible[key])for(const v of d[key])peak=Math.max(peak,v||0);for(const key of ['plus','minus']){if(!this.visible[key])continue;ctx.strokeStyle=this.colors[key];ctx.globalAlpha=.65;ctx.beginPath();let connected=false;for(let i=0;i<d[key].length;i++){if(d[key][i]===null){connected=false;continue}const x=(d.bin_start[i]+d.bin_end[i])/2/this.length*w,y=h-8-(this.scale(d[key][i])/this.scale(peak))*(h-18);if(connected)ctx.lineTo(x,y);else ctx.moveTo(x,y);connected=true}ctx.stroke()}ctx.globalAlpha=1;const a=this.range[0]/this.length*w,b=this.range[1]/this.length*w;ctx.fillStyle='#64748b18';ctx.fillRect(0,0,a,h);ctx.fillRect(b,0,w-b,h);ctx.fillStyle='#2563eb12';ctx.fillRect(a,0,Math.max(2,b-a),h);ctx.strokeStyle='#5685da';ctx.strokeRect(a+.5,.5,Math.max(1,b-a-1),h-1);ctx.fillStyle='#2563eb';ctx.fillRect(a,0,3,h);ctx.fillRect(Math.max(a,b-3),0,3,h)}
  showTooltip(event){if(!this.data)return;const d=this.data,pos=Math.min(this.range[1]-1,Math.floor(this.range[0]+this.hover*(this.range[1]-this.range[0])));let lo=0,hi=d.bin_end.length;while(lo<hi){const m=(lo+hi)>>1;if(d.bin_end[m]<=pos)lo=m+1;else hi=m}if(lo===d.bin_end.length||pos<d.bin_start[lo]){this.tooltip.hidden=true;return}const i=lo,range=d.bin_end[i]-d.bin_start[i]===1?`Position ${pos.toLocaleString()}`:`${d.bin_start[i].toLocaleString()}–${(d.bin_end[i]-1).toLocaleString()} · peak counts`;
    this.tooltip.replaceChildren();const title=document.createElement('strong');title.textContent=range+(d.dna?' · '+d.dna[pos-d.start]:'');this.tooltip.append(title);for(const key of ['plus','minus'])if(this.visible[key]){const line=document.createElement('div');line.style.color=this.colors[key];line.textContent=`csRNA ${key==='plus'?'+':'−'}: ${d[key][i]===null?'null (excluded)':d[key][i].toLocaleString()}`;this.tooltip.append(line)}const mask=document.createElement('div');mask.style.color='#72839a';mask.textContent=`${d.eligible[i].toLocaleString()} / ${(d.bin_end[i]-d.bin_start[i]).toLocaleString()} eligible bases`;this.tooltip.append(mask);this.tooltip.hidden=false;const r=this.canvas.getBoundingClientRect();this.tooltip.style.left=Math.max(0,Math.min(r.width-this.tooltip.offsetWidth-8,event.clientX-r.left+16))+'px';this.tooltip.style.top=Math.max(8,event.clientY-r.top-this.tooltip.offsetHeight-15)+'px';
  }
  save(format){this.hover=null;this.draw();let url;if(format==='png')url=this.canvas.toDataURL('image/png');else {const w=this.canvas.clientWidth,h=this.canvas.clientHeight,g=this.geometry(w,h);let content=`<rect width="100%" height="100%" fill="white"/><text x="72" y="18" font-family="sans-serif" font-size="12">csRNA counts${this.logScale?' · log10(count + 1)':''} · blue: + / orange: −</text>`;for(const letter of this.nucleotides(g)){content+=`<text x="${letter.x}" y="35" text-anchor="middle" font-family="monospace" font-size="12" font-weight="bold" fill="${letter.color}">${letter.base}</text>`}for(const value of this.ticks(g.peak)){const y=g.y(value);content+=`<path d="M72 ${y}H${g.right}" stroke="#edf1f6"/><text x="60" y="${y+4}" text-anchor="end" font-size="11">${value.toLocaleString(undefined,{maximumFractionDigits:1})}</text>`}for(const key of ['plus','minus']){if(!this.visible[key]||!this.data)continue;let path='',connected=false;const d=this.data;for(let i=0;i<d[key].length;i++){if(d.bin_end[i]<=this.range[0]||d.bin_start[i]>=this.range[1]||d[key][i]===null){connected=false;continue}const x=g.x((d.bin_start[i]+d.bin_end[i]-1)/2),y=g.y(d[key][i]);path+=`${connected?'L':'M'}${x},${y} `;connected=true}content+=`<path d="${path}" fill="none" stroke="${this.colors[key]}" stroke-width="1.6"/>`}content+=`<text x="72" y="${h-10}" font-family="sans-serif" font-size="12">Reference coordinates: ${this.range[0]}–${this.range[1]} (end excluded)</text>`;url=URL.createObjectURL(new Blob([`<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${content}</svg>`],{type:'image/svg+xml'}))}const a=document.createElement('a');a.href=url;a.download=`csrna-${this.range[0]}-${this.range[1]}.${format}`;a.click();if(format==='svg')setTimeout(()=>URL.revokeObjectURL(url),1000)}
}
