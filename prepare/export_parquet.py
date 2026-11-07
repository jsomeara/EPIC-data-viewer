"""Export seven coordinate-native Parquet files: one row per reference base."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time
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

def export_one(root,out,meta,batch_size):
    path=out/(meta['id']+'.parquet');record=out/(meta['id']+'.json')
    if path.exists() and record.exists():return json.loads(record.read_text())
    partial=path.with_suffix('.partial');dna_index=json.loads((root/(meta['assembly']+'.dna.json')).read_text())
    writer=None;rows=0;eligible_count=0;sums=[0,0];nonzero=[0,0];checks=[];last=time.monotonic()
    print('Exporting',meta['id'],'bases',sum(c['length'] for c in meta['contigs']),flush=True)
    try:
        with (root/(meta['assembly']+'.dna')).open('rb') as reference:
            for c in meta['contigs']:
                folder=root/meta['id']/c['id'];mask=np.fromfile(folder/'eligible.bin',dtype='<u8').reshape(-1,2)
                size=(folder/'signal.bin').stat().st_size
                signal=np.memmap(folder/'signal.bin',mode='r',dtype=SIGNAL) if size else np.empty(0,dtype=SIGNAL)
                for j,key in enumerate(('plus','minus')):
                    sums[j]+=int(signal[key].sum());nonzero[j]+=int(np.count_nonzero(signal[key]))
                eligible_count+=int((mask[:,1]-mask[:,0]).sum())
                reference.seek(dna_index[c['name']][0])
                for start in range(0,c['length'],batch_size):
                    end=min(start+batch_size,c['length']);dna=reference.read(end-start)
                    if len(dna)!=end-start:raise ValueError('Truncated DNA')
                    table=batch(meta,c,mask,signal,dna,start,end)
                    if writer is None:
                        writer=pq.ParquetWriter(partial,table.schema,compression='zstd',compression_level=3,
                                               use_dictionary=['assembly','organism','split','sequence_id','base'],
                                               column_encoding={'position':'DELTA_BINARY_PACKED'},write_page_checksum=True)
                    writer.write_table(table,row_group_size=len(table));rows+=len(table)
                    # Hash original logical row-group buffers for a full read-back comparison.
                    checks.append(table_digest(table))
                    if time.monotonic()-last>15:
                        print(meta['id'],rows,'bases written',flush=True);last=time.monotonic()
        if writer:writer.close();writer=None
        file=pq.ParquetFile(partial,page_checksum_verification=True)
        assert file.metadata.num_rows==rows
        print('Verifying',meta['id'],file.num_row_groups,'row groups',flush=True)
        for i,expected in enumerate(checks):
            if table_digest(file.read_row_group(i))!=expected:raise ValueError(f'Read-back mismatch: group {i}')
        del file
        partial.replace(path)
        result={'file':path.name,'dataset':meta['id'],'rows':rows,'eligible_bases':eligible_count,
                'count_sum_plus':sums[0],'count_sum_minus':sums[1],'nonzero_plus':nonzero[0],'nonzero_minus':nonzero[1],
                'bytes':path.stat().st_size,'sha256':digest(path),'row_groups':len(checks),'all_row_groups_verified':True}
        record.write_text(json.dumps(result,indent=2));print('Complete',path.name,result['bytes'],'bytes',flush=True)
        return result
    finally:
        if writer:writer.close()

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
    parser.add_argument('--batch-size',type=int,default=1_000_000);parser.add_argument('--workers',type=int,default=1);args=parser.parse_args()
    if args.batch_size<1:parser.error('batch size must be positive')
    args.output.mkdir(parents=True,exist_ok=True);catalog=json.loads((args.fixed/'catalog.json').read_text())
    jobs=[meta for meta in catalog if not args.dataset or meta['id'] in args.dataset]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(export_one,args.fixed,args.output,meta,args.batch_size) for meta in jobs]
        for future in futures:future.result()
    manifest={'format':'EPIC coordinate-native per-base v1','files':[json.loads(p.read_text()) for p in sorted(args.output.glob('*.json')) if p.name!='manifest.json']}
    temporary=args.output/'manifest.json.partial'
    temporary.write_text(json.dumps(manifest,indent=2))
    temporary.replace(args.output/'manifest.json')
if __name__=='__main__':main()
