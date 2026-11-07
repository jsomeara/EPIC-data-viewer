"""Export seven coordinate-native Parquet files: one row per reference base."""
import argparse
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time
import os
import threading
from functools import lru_cache
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SIGNAL=np.dtype([('position','<u8'),('plus','<u8'),('minus','<u8')])
SCHEMA=pa.schema([
    pa.field('assembly',pa.string(),nullable=False),pa.field('organism',pa.string(),nullable=False),
    pa.field('split',pa.string(),nullable=False),pa.field('sequence_id',pa.string(),nullable=False),
    pa.field('position',pa.uint64(),nullable=False),pa.field('base',pa.string(),nullable=False),
    pa.field('eligible',pa.bool_(),nullable=False),pa.field('csrna_plus',pa.uint64()),pa.field('csrna_minus',pa.uint64()),
    pa.field('label_plus',pa.uint8()),pa.field('label_minus',pa.uint8()),
],metadata={b'format':b'EPIC coordinate-native per-base v1',b'coordinates':b'zero-based; both strands follow reference orientation',
            b'nulls':b'excluded bases have null counts/labels; eligible measured zeros remain zero',
            b'rows':b'one row per reference base on each contig in the dataset whitelist; no artificial segments',
            b'targets':b'original supplied pooled csRNA counts; label = count > 0'})

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def constant(value,n):
    return pa.DictionaryArray.from_arrays(pa.array(np.zeros(n,dtype=np.int8)),pa.array([value]))

def batch(meta,c,mask,signal,dna,start,end):
    positions=np.arange(start,end,dtype=np.uint64);n=end-start
    idx=np.searchsorted(mask[:,1],positions,side='right');valid=idx<len(mask)
    eligible=valid&(mask[np.minimum(idx,len(mask)-1),0]<=positions)
    a,b=np.searchsorted(signal['position'],[start,end]);local=signal[a:b]
    plus=np.zeros(n,dtype=np.uint64);minus=plus.copy()
    offsets=(local['position']-start).astype(np.int64)
    plus[offsets]=local['plus'];minus[offsets]=local['minus']
    if np.any(~eligible[offsets]):raise ValueError('Signal outside whitelist')
    letters=np.frombuffer(dna,dtype=np.uint8)
    # Dictionary construction avoids allocating one Python object per nucleotide.
    bases=pa.DictionaryArray.from_arrays(pa.array(letters),pa.array([chr(i) for i in range(256)]))
    columns=[constant(meta['assembly'],n),constant(meta['organism'],n),constant(meta['split'],n),constant(c['name'],n),
             pa.array(positions),bases,pa.array(eligible),pa.array(plus,mask=~eligible),pa.array(minus,mask=~eligible),
             pa.array((plus>0).astype(np.uint8),mask=~eligible),pa.array((minus>0).astype(np.uint8),mask=~eligible)]
    # Keep dictionary arrays until the writer rather than materializing repeated strings.
    schema=pa.schema([pa.field(field.name,col.type,nullable=field.nullable) for field,col in zip(SCHEMA,columns)],metadata=SCHEMA.metadata)
    return pa.Table.from_arrays(columns,schema=schema)

def export_one(root,out,meta,batch_size,verify_workers=4):
    path=out/(meta['id']+'.parquet');record=out/(meta['id']+'.json')
    if path.exists() and record.exists():return json.loads(record.read_text())
    partial=path.with_suffix('.partial');dna_index=json.loads((root/(meta['assembly']+'.dna.json')).read_text())
    groups=[(c,start,min(start+batch_size,c['length'])) for c in meta['contigs'] for start in range(0,c['length'],batch_size)]
    rows=sum(c['length'] for c in meta['contigs']);complete=False
    if partial.exists():
        try:
            existing=pq.ParquetFile(partial)
            complete=existing.metadata.num_rows==rows and existing.num_row_groups==len(groups)
            del existing
        except pa.ArrowInvalid:pass
    @lru_cache(maxsize=32)
    def arrays(contig_id):
        folder=root/meta['id']/contig_id
        mask=np.fromfile(folder/'eligible.bin',dtype='<u8').reshape(-1,2)
        signal=np.memmap(folder/'signal.bin',mode='r',dtype=SIGNAL) if (folder/'signal.bin').stat().st_size else np.empty(0,dtype=SIGNAL)
        return mask,signal
    writer=None;last=time.monotonic();written=0
    print(('Resuming verification' if complete else 'Exporting'),meta['id'],'bases',rows,flush=True)
    with (root/(meta['assembly']+'.dna')).open('rb') as reference:
        try:
            if not complete:
                for c,start,end in groups:
                    mask,signal=arrays(c['id'])
                    dna=os.pread(reference.fileno(),end-start,dna_index[c['name']][0]+start)
                    if len(dna)!=end-start:raise ValueError('Truncated DNA')
                    table=batch(meta,c,mask,signal,dna,start,end)
                    if writer is None:
                        writer=pq.ParquetWriter(partial,table.schema,compression='zstd',compression_level=3,
                                               use_dictionary=['assembly','organism','split','sequence_id','base'],
                                               column_encoding={'position':'DELTA_BINARY_PACKED'},write_page_checksum=True)
                    writer.write_table(table,row_group_size=len(table));written+=len(table)
                    if time.monotonic()-last>15:
                        print(meta['id'],written,'bases written',flush=True);last=time.monotonic()
                if writer:writer.close();writer=None
            # Independent readers avoid sharing mutable reader state between threads.
            # Reconstruct every original row group from coordinate indexes, then compare
            # all logical columns; batching and dictionary order do not alter meaning.
            check=pq.ParquetFile(partial)
            if check.metadata.num_rows!=rows or check.num_row_groups!=len(groups):raise ValueError('Wrong Parquet row/group count')
            del check
            local=threading.local()
            def verify(i):
                if not hasattr(local,'file'):local.file=pq.ParquetFile(partial,page_checksum_verification=True)
                c,start,end=groups[i];mask,signal=arrays(c['id'])
                dna=os.pread(reference.fileno(),end-start,dna_index[c['name']][0]+start)
                expected=batch(meta,c,mask,signal,dna,start,end)
                actual=local.file.read_row_group(i,use_threads=False)
                if table_digest(actual)!=table_digest(expected):raise ValueError(f'Read-back mismatch: {meta["id"]} group {i}')
            print('Verifying',meta['id'],len(groups),'row groups with',verify_workers,'threads',flush=True)
            with ThreadPoolExecutor(max_workers=verify_workers) as pool:
                for _ in pool.map(verify,range(len(groups))):pass
        finally:
            if writer:writer.close()
    partial.replace(path)
    sums=[0,0];nonzero=[0,0];eligible_count=0
    for c in meta['contigs']:
        mask,signal=arrays(c['id']);eligible_count+=int((mask[:,1]-mask[:,0]).sum())
        for j,key in enumerate(('plus','minus')):
            sums[j]+=int(signal[key].sum());nonzero[j]+=int(np.count_nonzero(signal[key]))
    result={'file':path.name,'dataset':meta['id'],'rows':rows,'eligible_bases':eligible_count,
            'count_sum_plus':sums[0],'count_sum_minus':sums[1],'nonzero_plus':nonzero[0],'nonzero_minus':nonzero[1],
            'bytes':path.stat().st_size,'sha256':digest(path),'row_groups':len(groups),'all_row_groups_verified':True}
    record.write_text(json.dumps(result,indent=2));print('Complete',path.name,result['bytes'],'bytes',flush=True)
    return result

def table_digest(table):
    import pyarrow.compute as pc
    h=hashlib.sha256();h.update(str(len(table)).encode())
    for name in table.column_names:
        col=table[name].combine_chunks()
        if name in ('assembly','organism','split','sequence_id'):
            unique=pc.unique(col).to_pylist()
            if len(unique)!=1:raise ValueError('A row group crosses a contig or dataset')
            h.update(unique[0].encode())
        elif name=='base':
            if not pa.types.is_dictionary(col.type):col=pc.dictionary_encode(col)
            lookup=np.array([ord(v) for v in col.dictionary.to_pylist()],dtype=np.uint8)
            h.update(lookup[col.indices.to_numpy(zero_copy_only=False)].tobytes())
        elif name in ('csrna_plus','csrna_minus','label_plus','label_minus'):
            h.update(pc.is_valid(col).to_numpy(zero_copy_only=False).tobytes())
            h.update(pc.fill_null(col,0).to_numpy(zero_copy_only=False).tobytes())
        else:h.update(col.to_numpy(zero_copy_only=False).tobytes())
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--fixed',type=Path,default=Path('data_fixed'))
    parser.add_argument('--output',type=Path,default=Path('parquet'));parser.add_argument('--dataset',action='append')
    parser.add_argument('--batch-size',type=int,default=1_000_000);parser.add_argument('--workers',type=int,default=1);parser.add_argument('--verify-workers',type=int,default=4);args=parser.parse_args()
    if args.batch_size<1:parser.error('batch size must be positive')
    args.output.mkdir(parents=True,exist_ok=True);catalog=json.loads((args.fixed/'catalog.json').read_text())
    jobs=[meta for meta in catalog if not args.dataset or meta['id'] in args.dataset]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(export_one,args.fixed,args.output,meta,args.batch_size,args.verify_workers) for meta in jobs]
        for future in futures:future.result()
    manifest={'format':'EPIC coordinate-native per-base v1','files':[json.loads(p.read_text()) for p in sorted(args.output.glob('*.json')) if p.name!='manifest.json']}
    temporary=args.output/'manifest.json.partial'
    temporary.write_text(json.dumps(manifest,indent=2))
    temporary.replace(args.output/'manifest.json')
if __name__=='__main__':main()
