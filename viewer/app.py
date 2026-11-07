"""Read-only coordinate browser for data_fixed. No segment boundaries."""
import argparse
import gzip
import json
import logging
import math
import re
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
import numpy as np

ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/'data_fixed'
DOWNLOADS=ROOT.parent/'parquet'
SIGNAL=np.dtype([('position','<u8'),('plus','<u8'),('minus','<u8')])

@lru_cache(maxsize=1)
def catalog():
    return {m['id']:m for m in json.loads((DATA/'catalog.json').read_text())}

@lru_cache(maxsize=6)
def reference_index(assembly):
    return json.loads((DATA/(assembly+'.dna.json')).read_text())

@lru_cache(maxsize=32)
def track(dataset,contig):
    meta=catalog()[dataset]
    c=next((c for c in meta['contigs'] if c['id']==contig),None)
    if c is None:raise ValueError('Unknown chromosome/contig')
    folder=DATA/dataset/contig
    intervals=np.fromfile(folder/'eligible.bin',dtype='<u8').reshape(-1,2)
    signal=read_array(folder/'signal.bin')
    dna_index=reference_index(meta['assembly'])
    dna_offset,dna_length=dna_index[c['name']]
    assert dna_length==c['length']
    return c,folder,intervals,signal,dna_offset,DATA/(meta['assembly']+'.dna')

def read_array(path):
    return np.memmap(path,dtype=SIGNAL,mode='r') if path.stat().st_size else np.empty(0,dtype=SIGNAL)

@lru_cache(maxsize=64)
def peaks(dataset,contig,width):
    return read_array(DATA/dataset/contig/f'peaks-{width}.bin')

def bounds(dataset,contig,start,end):
    c,*_=track(dataset,contig)
    if not 0<=start<end<=c['length']:raise ValueError('Use 0 ≤ start < end ≤ contig length')

def eligibility(intervals,positions):
    indices=np.searchsorted(intervals[:,1],positions,side='right')
    clipped=np.minimum(indices,len(intervals)-1)
    return (indices<len(intervals)) & (intervals[clipped,0]<=positions)

def coverage(intervals,edges):
    # Integral of the eligibility mask at each edge, without expanding per-base data.
    lengths=intervals[:,1]-intervals[:,0]
    prefix=np.r_[np.uint64(0),np.cumsum(lengths)]
    indices=np.searchsorted(intervals[:,1],edges,side='right')
    clipped=np.minimum(indices,len(intervals)-1)
    partial=np.maximum(0,edges.astype(np.int64)-intervals[clipped,0].astype(np.int64))
    partial=np.minimum(partial,lengths[clipped])
    cumulative=prefix[indices]+np.where(indices<len(intervals),partial,0)
    return np.diff(cumulative).astype(np.int64)

def slice_signal(signal,start,end):
    a,b=np.searchsorted(signal['position'],[start,end])
    return signal[a:b]

def view(dataset,contig,start,end,pixels=1200):
    bounds(dataset,contig,start,end)
    c,folder,intervals,signal,*_=track(dataset,contig)
    points=max(200,min(2000,pixels))
    desired=math.ceil((end-start)/points)
    width=1 if desired==1 else desired
    if desired>=256:
        width=256
        while width<desired:width*=4
        # Small tracks can end their pyramid early; the last level is sufficient.
        source_width=width
        while not (folder/f'peaks-{source_width}.bin').exists() and source_width>256:source_width//=4
        if not (folder/f'peaks-{source_width}.bin').exists():source=signal
        else:source=peaks(dataset,contig,source_width)
    else:source=signal
    left=(start//width)*width
    edges=np.arange(left,end+width,width,dtype=np.int64)
    n=len(edges)-1
    plus=np.zeros(n,dtype=np.uint64);minus=plus.copy()
    data=slice_signal(source,left,int(edges[-1]))
    # Edge bins must be computed from exact data, not peaks outside the viewport.
    indices=((data['position'].astype(np.int64)-left)//width)
    for dest,key in ((plus,'plus'),(minus,'minus')):
        np.maximum.at(dest,indices,data[key])
    for i in {0,n-1}:
        exact=slice_signal(signal,max(start,int(edges[i])),min(end,int(edges[i+1])))
        plus[i]=exact['plus'].max(initial=0);minus[i]=exact['minus'].max(initial=0)
    clipped=np.clip(edges,start,end)
    eligible=coverage(intervals,clipped)
    def channel(values):return [int(v) if eligible[i] else None for i,v in enumerate(values)]
    dna=None
    if width==1:
        *_,dna_offset,dna_path=track(dataset,contig)
        with dna_path.open('rb') as stream:
            stream.seek(dna_offset+start);dna=stream.read(end-start).decode()
    return {'start':start,'end':end,'length':c['length'],'bin_size':width,'dna':dna,
            'bin_start':clipped[:-1].tolist(),'bin_end':clipped[1:].tolist(),
            'plus':channel(plus),'minus':channel(minus),'eligible':eligible.tolist(),
            'has_signal':bool(np.any(plus)|np.any(minus))}

def base_rows(dataset,contig,start,end,cursor=None,nonzero=False,limit=100):
    bounds(dataset,contig,start,end)
    c,_,intervals,signal,dna_offset,dna_path=track(dataset,contig)
    cursor=start if cursor is None else cursor
    if not start<=cursor<=end:raise ValueError('Invalid page cursor')
    if nonzero:
        values=slice_signal(signal,cursor,end)[:limit+1]
        positions=values['position'][:limit].astype(np.int64)
        more=len(values)>limit
    else:
        positions=np.arange(cursor,min(end,cursor+limit),dtype=np.int64)
        more=cursor+limit<end
    if not len(positions):return {'rows':[],'next':None}
    a,b=np.searchsorted(signal['position'],[positions[0],positions[-1]+1])
    local=signal[a:b]
    indices=np.searchsorted(local['position'],positions)
    counts=np.zeros((len(positions),2),dtype=np.uint64)
    valid=indices<len(local)
    matches=np.flatnonzero(valid)
    matches=matches[local['position'][indices[matches]]==positions[matches]]
    counts[matches,0]=local['plus'][indices[matches]];counts[matches,1]=local['minus'][indices[matches]]
    mask=eligibility(intervals,positions)
    result=[]
    with dna_path.open('rb') as dna:
        # Non-zero pages may span millions of bases; read only the requested letters.
        if nonzero:
            letters=[]
            for p in positions:dna.seek(dna_offset+int(p));letters.append(dna.read(1).decode())
        else:
            dna.seek(dna_offset+int(positions[0]));letters=dna.read(len(positions)).decode()
    for i,p in enumerate(positions):
        plus,minus=(int(v) for v in counts[i])
        result.append([int(p),letters[i],bool(mask[i]),plus if mask[i] else None,minus if mask[i] else None,
                       int(plus>0) if mask[i] else None,int(minus>0) if mask[i] else None])
    return {'rows':result,'next':int(positions[-1])+1 if more else None}

def downloads():
    manifest=DOWNLOADS/'manifest.json'
    if not manifest.exists():return []
    result=[]
    for record in json.loads(manifest.read_text())['files']:
        name=record['file']
        if not re.fullmatch(r'[a-z0-9_]+\.parquet',name):continue
        path=DOWNLOADS/name
        if path.is_file() and path.stat().st_size==record['bytes']:
            meta=catalog()[record['dataset']]
            result.append({**record,'organism':meta['organism'],'split':meta['split'],'url':'downloads/'+name})
    return result

class Handler(BaseHTTPRequestHandler):
    def do_HEAD(self):
        self.do_GET(head_only=True)

    def serve_download(self,path,head_only):
        name=path.removeprefix('/downloads/')
        allowed={record['file'] for record in downloads()}|{'manifest.json'}
        if name not in allowed or not (DOWNLOADS/name).is_file():self.send_error(404);return
        file=DOWNLOADS/name;size=file.stat().st_size;start=0;end=size-1;status=200
        if self.headers.get('Range'):
            match=re.fullmatch(r'bytes=(\d*)-(\d*)',self.headers['Range'])
            if not match or not any(match.groups()):self.send_error(416);return
            if match[1]:start=int(match[1]);end=min(size-1,int(match[2])) if match[2] else size-1
            else:start=max(0,size-int(match[2]))
            if start> end or start>=size:
                self.send_response(416);self.send_header('Content-Range',f'bytes */{size}');self.end_headers();return
            status=206
        self.send_response(status)
        self.send_header('Content-Type','application/vnd.apache.parquet' if name.endswith('.parquet') else 'application/json')
        self.send_header('Content-Disposition',f'attachment; filename="{name}"')
        self.send_header('Content-Length',str(end-start+1));self.send_header('Accept-Ranges','bytes')
        if status==206:self.send_header('Content-Range',f'bytes {start}-{end}/{size}')
        self.end_headers()
        if head_only:return
        with file.open('rb') as stream:
            stream.seek(start);remaining=end-start+1
            while remaining:
                block=stream.read(min(256*1024,remaining))
                if not block:break
                self.wfile.write(block);remaining-=len(block)

    def do_GET(self,head_only=False):
        req=urlparse(self.path)
        try:
            if req.path.startswith('/downloads/'):
                self.serve_download(req.path,head_only);return
            assets={'/':'index.html','/plot.js':'plot.js','/viewer.js':'viewer.js','/style.css':'style.css'}
            if req.path in assets:
                file=ROOT/assets[req.path];body=file.read_bytes()
                kind={'.html':'text/html','.js':'text/javascript','.css':'text/css'}[file.suffix]+'; charset=utf-8'
            elif req.path.startswith('/api/'):
                q={k:v[0] for k,v in parse_qs(req.query).items()}
                if req.path=='/api/downloads':result=downloads()
                elif req.path=='/api/files':
                    result=[{k:v for k,v in m.items() if k!='contigs'} for m in catalog().values()]
                elif req.path=='/api/contigs':result=catalog()[q['dataset']]['contigs']
                elif req.path in ('/api/view','/api/bases'):
                    d,c=q['dataset'],q['contig'];start,end=int(q['start']),int(q['end'])
                    if req.path=='/api/view':result=view(d,c,start,end,int(q.get('pixels','1200')))
                    else:result=base_rows(d,c,start,end,int(q['cursor']) if 'cursor' in q else None,q.get('nonzero')=='1')
                else:self.send_error(404);return
                body=json.dumps(result,separators=(',',':')).encode();kind='application/json'
            else:self.send_error(404);return
            self.send_response(200)
            if 'gzip' in self.headers.get('Accept-Encoding','') and len(body)>2000:
                body=gzip.compress(body,compresslevel=1);self.send_header('Content-Encoding','gzip')
                self.send_header('Vary','Accept-Encoding')
            self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(body)))
            self.send_header('Cache-Control','no-cache');self.end_headers()
            if not head_only:self.wfile.write(body)
        except (ValueError,KeyError) as error:self.send_error(400,'Invalid request',explain=str(error))
        except (BrokenPipeError,ConnectionResetError):pass
        except Exception:
            logging.exception('Viewer request failed');self.send_error(500,'Unable to read data_fixed; check server log')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8767)
    parser.add_argument('--host',default='127.0.0.1');args=parser.parse_args()
    if not (DATA/'catalog.json').exists():parser.error('data_fixed/catalog.json is missing. Follow the repository preparation instructions before starting the viewer.')
    server=ThreadingHTTPServer((args.host,args.port),Handler)
    print(f'Coordinate viewer: http://{args.host}:{args.port}',flush=True);server.serve_forever()
