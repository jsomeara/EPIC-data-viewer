"""Build segment-free viewer indexes directly from the extracted EPIC Zenodo data."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import numpy as np

SIGNAL=np.dtype([('position','<u8'),('plus','<u8'),('minus','<u8')])
DATASETS={
 'ASM119413v2':'california_two_spot_octopus', 'Ofas_2.0':'milkweed_bug',
 'cgigas_uk_roslin_v1':'pacific_oyster','ilPloInte3.2':'indianmeal_moth',
 'sSquSuc5.hap2':'puget_sound_dogfish','jaNemVect1.1':'starlet_sea_anemone_toy'}

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def lines(path):
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.strip():yield line.split()

def reference(raw,output,assembly):
    source=raw/assembly/'genome';file=output/(assembly+'.dna');indexfile=output/(assembly+'.dna.json')
    if file.exists() and indexfile.exists():return json.loads(indexfile.read_text())
    expected={v[0]:int(v[1]) for v in (line.split() for line in (source/'genome.fa.fai').read_text().splitlines())}
    index={};name=None;offset=0;length=0;pending=b''
    with gzip.open(source/'genome.fa.gz','rb') as genome,file.with_suffix('.partial').open('wb') as dest:
        def write(piece):
            nonlocal length
            dna=piece.replace(b'\n',b'').replace(b'\r',b'').upper()
            if name is None and dna:raise ValueError('DNA before FASTA header')
            dest.write(dna);length+=len(dna)
        while True:
            block=genome.read(8*1024*1024)
            if not block:break
            pending+=block;edge=pending.rfind(b'\n')+1
            if not edge:continue
            complete,pending=pending[:edge],pending[edge:];cursor=0
            for match in re.finditer(rb'^>([^\s]+)[^\n]*\n',complete,re.MULTILINE):
                write(complete[cursor:match.start()])
                if name is not None:
                    if name in index or expected[name]!=length:raise ValueError('FASTA index mismatch: '+name)
                    index[name]=[offset,length]
                name=match[1].decode();offset=dest.tell();length=0;cursor=match.end()
            write(complete[cursor:])
        write(pending)
        if name is not None:
            if name in index or expected[name]!=length:raise ValueError('FASTA index mismatch: '+name)
            index[name]=[offset,length]
    if set(index)!=set(expected):raise ValueError('FASTA records do not match .fai')
    file.with_suffix('.partial').replace(file);indexfile.write_text(json.dumps(index))
    return index

def vector_nonzero(path):
    """Parse dense vectors in bounded chunks, yielding only positive entries."""
    offset=0;pending=b''
    with gzip.open(path,'rb') as stream:
        while True:
            block=stream.read(8*1024*1024)
            if not block:
                if not pending:break
                block=b'\n'
            pending+=block;edge=pending.rfind(b'\n')+1
            if not edge:continue
            complete,pending=pending[:edge],pending[edge:]
            values=np.fromstring(complete.decode('ascii'),dtype=np.float64,sep='\n')
            if len(values)!=complete.count(b'\n'):raise ValueError('Malformed target vector')
            if not np.all(np.isfinite(values)) or np.any(values<0) or np.any(values>2**53) or np.any(values!=np.floor(values)):
                raise ValueError('Targets must be nonnegative exactly representable integer counts')
            nz=np.flatnonzero(values)
            yield offset+nz,values[nz].astype(np.uint64),len(values)
            offset+=len(values)

def pyramid(source,folder):
    width=256
    while len(source):
        bins=source['position']//width;first=np.r_[0,np.flatnonzero(bins[1:]!=bins[:-1])+1]
        reduced=np.empty(len(first),dtype=SIGNAL);reduced['position']=bins[first]*width
        for key in ('plus','minus'):reduced[key]=np.maximum.reduceat(source[key],first)
        reduced.tofile(folder/f'peaks-{width}.bin')
        if len(reduced)<=1:break
        source=reduced;width*=4

def build(raw,output,assembly,split):
    name=DATASETS[assembly]+'_'+split;folder=output/name;meta_path=folder/'meta.json'
    if meta_path.exists():return json.loads(meta_path.read_text())
    folder.mkdir(parents=True,exist_ok=True);source=raw/assembly
    index=reference(raw,output,assembly)
    organism=next(line.split(':',1)[1].strip().split(' (')[0] for line in (source/'genome/assembly_report.txt').read_text().splitlines() if line.startswith('# Organism name:'))
    whitelist=[(v[0],int(v[1]),int(v[2])) for v in lines(source/f'for_participants/whitelist.{split}.bed.gz')]
    if not whitelist:raise ValueError('Empty whitelist')
    by_seq={}
    for seq,start,end in whitelist:
        if not 0<=start<end<=index[seq][1]:raise ValueError('Invalid whitelist interval')
        by_seq.setdefault(seq,[]).append((start,end))
    contigs=[];by_name={}
    for i,(seq,spans) in enumerate(by_seq.items()):
        dest=folder/str(i);dest.mkdir(exist_ok=True);mask=np.asarray(spans,dtype='<u8')
        if np.any(mask[1:,0]<mask[:-1,1]):raise ValueError('Overlapping/unsorted whitelist')
        mask.tofile(dest/'eligible.bin');(dest/'signal.bin').write_bytes(b'')
        c={'id':str(i),'name':seq,'length':index[seq][1],'eligible_bases':int((mask[:,1]-mask[:,0]).sum())}
        contigs.append(c);by_name[seq]=c
    # Validate the exact vector order against the provided template, not genomic sorting.
    template=iter(lines(source/f'for_participants/template.{split}.bed.gz'))
    for strand in ('+','-'):
        for seq,start,end in whitelist:
            v=next(template,None)
            if v is None or (v[0],int(v[1]),int(v[2]),v[5])!=(seq,start,end,strand):
                raise ValueError('Template order differs from whitelist-plus then whitelist-minus')
    if next(template,None) is not None:raise ValueError('Unexpected extra template entries')
    starts=np.asarray([v[1] for v in whitelist],dtype=np.uint64)
    spans=np.asarray([v[2]-v[1] for v in whitelist],dtype=np.uint64)
    ends=np.cumsum(spans);begins=np.r_[np.uint64(0),ends[:-1]];half=int(ends[-1])
    vector=source/(f'for_participants/csRNA-{split}.txt.gz' if split=='train' else 'ground_truth/csRNA-test.txt.gz')
    total=0;source_sum=0;source_nonzero=0
    print('Reading pooled targets',assembly,split,flush=True)
    for indices,counts,length in vector_nonzero(vector):
        total+=length;source_sum+=int(counts.sum());source_nonzero+=len(counts)
        if len(indices) and int(indices[-1])>=2*half:raise ValueError('Target vector is too long')
        strand=indices>=half;within=indices%half;interval=np.searchsorted(ends,within,side='right')
        position=starts[interval]+within-begins[interval]
        groups={}
        for j,k in enumerate(interval):groups.setdefault(whitelist[int(k)][0],[]).append(j)
        for seq,selection in groups.items():
            selection=np.asarray(selection);records=np.zeros(len(selection),dtype=SIGNAL)
            records['position']=position[selection]
            records['plus']=np.where(strand[selection],0,counts[selection])
            records['minus']=np.where(strand[selection],counts[selection],0)
            with (folder/by_name[seq]['id']/'signal.bin').open('ab') as out:records.tofile(out)
    if total!=2*half:raise ValueError(f'Wrong target vector length: {total} != {2*half}')
    final_sum=0;final_nonzero=0
    for c in contigs:
        dest=folder/c['id'];array=np.fromfile(dest/'signal.bin',dtype=SIGNAL)
        if len(array):
            array.sort(order='position');first=np.r_[0,np.flatnonzero(array['position'][1:]!=array['position'][:-1])+1]
            merged=np.empty(len(first),dtype=SIGNAL);merged['position']=array['position'][first]
            for key in ('plus','minus'):
                if np.any(np.add.reduceat((array[key]>0).astype(np.uint8),first)>1):raise ValueError('Duplicate base-and-strand target')
                merged[key]=np.add.reduceat(array[key],first)
                final_sum+=int(merged[key].sum());final_nonzero+=int(np.count_nonzero(merged[key]))
            merged.tofile(dest/'signal.bin');pyramid(merged,dest);c['nonzero_bases']=len(merged)
        else:c['nonzero_bases']=0
    if source_sum!=final_sum or source_nonzero!=final_nonzero:raise ValueError('Target totals changed during indexing')
    provenance=[]
    for path in (source/'genome/genome.fa.gz',source/'genome/genome.fa.fai',source/f'for_participants/whitelist.{split}.bed.gz',source/f'for_participants/template.{split}.bed.gz',vector):
        provenance.append({'path':str(path.relative_to(raw)),'bytes':path.stat().st_size,'sha256':sha256(path)})
    meta={'id':name,'assembly':assembly,'organism':organism,'split':split,'contigs':contigs,
          'source':str(vector.relative_to(raw)),'nonzero_bases':sum(c['nonzero_bases'] for c in contigs),
          'source_files':provenance,'source_nonzero_targets':source_nonzero,'source_count_sum':source_sum}
    meta_path.write_text(json.dumps(meta));print('Indexed',name,flush=True);return meta

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--raw',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path('data_fixed'));parser.add_argument('--assembly',choices=DATASETS,action='append');args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True);catalog=[]
    for assembly in args.assembly or DATASETS:
        for split in (('train','test') if assembly=='jaNemVect1.1' else ('train',)):
            catalog.append(build(args.raw,args.output,assembly,split))
    (args.output/'catalog.json').write_text(json.dumps(catalog))
if __name__=='__main__':main()
