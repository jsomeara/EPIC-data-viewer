const $=id=>document.getElementById(id);
const state={dataset:null,contigs:[],contig:null,start:0,end:1,cursor:0,pages:[],next:null};
const pending=new Map();let viewportTimer;
function error(message){$('error').textContent=message;$('error').hidden=!message}
async function request(slot,path,indicator){
  pending.get(slot)?.abort();const controller=new AbortController();pending.set(slot,controller);
  if(indicator){$(indicator).hidden=false;$(indicator).parentElement.setAttribute('aria-busy','true')}
  try{const response=await fetch(path,{signal:controller.signal});if(!response.ok)throw Error(`Could not read this region (HTTP ${response.status}). Check server.log.`);const data=await response.json();return pending.get(slot)===controller&&!controller.signal.aborted?data:null}
  catch(e){if(e.name!=='AbortError'){error(e.message)}return null}
  finally{if(pending.get(slot)===controller){pending.delete(slot);if(indicator){$(indicator).hidden=true;$(indicator).parentElement.setAttribute('aria-busy','false')}}}
}
function parameters(extra={}){return new URLSearchParams({dataset:state.dataset,contig:state.contig.id,start:state.start,end:state.end,...extra})}
const plot=new TrackPlot((start,end,debounce)=>{
  if(!state.contig)return;
  clearTimeout(viewportTimer);
  pending.get('view')?.abort();pending.get('table')?.abort();
  state.start=start;state.end=end;updateCoordinates();
  if(debounce)viewportTimer=setTimeout(refreshRegion,120);else refreshRegion();
});
function updateCoordinates(){
  $('start').value=state.start;$('end').value=state.end;
  $('locus').textContent=`${state.contig.name}:${state.start.toLocaleString()}–${state.end.toLocaleString()} · ${(state.end-state.start).toLocaleString()} bp`;
}
async function refreshRegion(){
  if(!state.contig)return;error('');state.cursor=state.start;state.pages=[];loadTable();
  const data=await request('view','api/view?'+parameters({pixels:Math.round($('plot').clientWidth)}),'chartLoading');
  if(!data)return;plot.setData(data);$('noSignal').hidden=data.has_signal;
  $('resolution').textContent=data.bin_size===1?'Exact per-base counts':`Peak counts per ${data.bin_size.toLocaleString()} bp bin · zoom in for exact bases`;
}
async function loadTable(){
  if(!state.contig)return;
  $('previousPage').disabled=true;$('nextPage').disabled=true;
  const data=await request('table','api/bases?'+parameters({cursor:state.cursor,nonzero:$('nonzero').checked?'1':'0'}),'tableLoading');
  if(!data)return;state.next=data.next;const fragment=document.createDocumentFragment();
  for(const values of data.rows){const row=document.createElement('tr');if(values[3]>0||values[4]>0)row.className='positive';for(let i=0;i<values.length;i++){const cell=document.createElement('td');cell.textContent=values[i]===null?'null':typeof values[i]==='number'?values[i].toLocaleString():values[i];if(values[i]===null)cell.className='null';row.append(cell)}fragment.append(row)}
  if(!data.rows.length){const row=document.createElement('tr'),cell=document.createElement('td');cell.colSpan=7;cell.textContent='No matching base pairs in this region.';cell.style.textAlign='center';row.append(cell);fragment.append(row)}
  $('bases').replaceChildren(fragment);$('bases').closest('.table-wrap').scrollTop=0;
  $('pageInfo').textContent=data.rows.length?`Page ${state.pages.length+1} · ${data.rows[0][0].toLocaleString()}–${data.rows.at(-1)[0].toLocaleString()} · ${data.rows.length} base pairs`:'No matches';
  $('previousPage').disabled=!state.pages.length;$('nextPage').disabled=data.next===null;
}
async function chooseContig(){
  clearTimeout(viewportTimer);
  pending.get('view')?.abort();pending.get('table')?.abort();
  state.contig=state.contigs.find(c=>c.id===$('contig').value);if(!state.contig)return;
  state.start=0;state.end=state.contig.length;plot.data=null;plot.overData=null;
  plot.setRange(state.start,state.end,state.contig.length);updateCoordinates();
  $('trackInfo').textContent=`${state.contig.length.toLocaleString()} bp · ${state.contig.eligible_bases.toLocaleString()} eligible · ${state.contig.nonzero_bases.toLocaleString()} non-zero`;
  $('start').max=state.contig.length-1;$('end').max=state.contig.length;
  refreshRegion();
  const data=await request('overview','api/view?'+parameters({start:0,end:state.contig.length,pixels:1200}),'overviewLoading');
  if(data)plot.setOverview(data);
}
async function chooseDataset(){
  clearTimeout(viewportTimer);
  state.dataset=$('dataset').value;state.contig=null;
  for(const slot of ['view','table','overview','nextSignal'])pending.get(slot)?.abort();
  const data=await request('contigs','api/contigs?'+new URLSearchParams({dataset:state.dataset}),'chartLoading');if(!data)return;
  state.contigs=data;$('contig').replaceChildren();
  for(const c of data){const option=document.createElement('option');option.value=c.id;option.textContent=`${c.name} · ${c.length.toLocaleString()} bp`;$('contig').append(option)}
  // Start on the largest chromosome/contig instead of a small scaffold.
  $('contig').value=data.reduce((largest,c)=>c.length>largest.length?c:largest,data[0]).id;
  chooseContig();
}
$('dataset').onchange=chooseDataset;$('contig').onchange=chooseContig;
$('regionForm').onsubmit=e=>{e.preventDefault();if(!state.contig)return;const a=Number($('start').value),b=Number($('end').value);if(!Number.isInteger(a)||!Number.isInteger(b)||a<0||b<=a||b>state.contig.length){error('Enter a valid interval: 0 ≤ start < end ≤ contig length.');return}plot.navigate(a,b)};
$('zoomIn').onclick=()=>plot.zoom(.5);$('zoomOut').onclick=()=>plot.zoom(2);
$('panLeft').onclick=()=>plot.shift(-.5);$('panRight').onclick=()=>plot.shift(.5);
$('reset').onclick=()=>plot.navigate(0,plot.length);
$('mode').onclick=()=>{plot.pan=!plot.pan;$('mode').textContent=plot.pan?'Drag: pan':'Drag: zoom';$('mode').setAttribute('aria-pressed',plot.pan);$('plot').style.cursor=plot.pan?'grab':'crosshair'};
for(const key of ['plus','minus'])$(key).onchange=()=>{plot.visible[key]=$(key).checked;plot.schedule()};
$('logScale').onchange=()=>{plot.logScale=$('logScale').checked;plot.schedule()};
$('export').onclick=()=>plot.save('png');$('exportSvg').onclick=()=>plot.save('svg');
$('nonzero').onchange=()=>{state.cursor=state.start;state.pages=[];loadTable()};
$('nextPage').onclick=()=>{if(state.next===null)return;state.pages.push(state.cursor);state.cursor=state.next;loadTable()};
$('previousPage').onclick=()=>{if(!state.pages.length)return;state.cursor=state.pages.pop();loadTable()};
$('nextSignal').onclick=async()=>{
  if(!state.contig)return;
  const length=state.contig.length;let start=state.end<length?state.end:state.start;
  let data=await request('nextSignal','api/bases?'+parameters({start,end:length,cursor:start,nonzero:1}),'chartLoading');
  if(!data)return;
  if(!data.rows.length&&start>0)data=await request('nextSignal','api/bases?'+parameters({start:0,end:start,cursor:0,nonzero:1}),'chartLoading');
  if(!data?.rows.length){error('No non-zero signal found on this contig.');return}
  const position=data.rows[0][0];plot.navigate(position-1000,position+1001);
};
request('files','api/files','chartLoading').then(data=>{if(!data)return;for(const file of data){const option=document.createElement('option');option.value=file.id;option.textContent=`${file.organism} · ${file.split}`;$('dataset').append(option)}chooseDataset()});

function loadDownloads(){request('downloads','api/downloads').then(files=>{
  if(!files)return;
  const fragment=document.createDocumentFragment();
  for(const file of files){
    const row=document.createElement('div');row.className='download-file';
    const link=document.createElement('a');link.href=file.url;link.download=file.file;
    link.textContent=`${file.organism} · ${file.split}`;link.title=file.file;
    const size=document.createElement('span');size.textContent=`${(file.bytes/1024/1024).toFixed(1)} MiB · ${file.rows.toLocaleString()} bases`;
    row.append(link,size);fragment.append(row);
  }
  if(!files.length){const text=document.createElement('p');text.className='muted';text.textContent='Parquet files are being prepared.';fragment.append(text)}
  $('downloads').replaceChildren(fragment);$('downloadManifest').hidden=!files.length;
});}
$('downloadMenu').ontoggle=()=>{if($('downloadMenu').open)loadDownloads()};
loadDownloads();
