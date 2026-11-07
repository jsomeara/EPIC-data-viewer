"""Synthetic raw-to-Parquet and download tests; no EPIC data required."""
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
import sys
from urllib.parse import urlencode
import pyarrow.parquet as pq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare.from_raw import build,vector_nonzero
from prepare.export_parquet import export_one
from viewer import app

class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.raw=self.root/'raw';self.fixed=self.root/'fixed';self.fixed.mkdir()
        source=self.raw/'jaNemVect1.1';(source/'genome').mkdir(parents=True);(source/'for_participants').mkdir()
        (source/'genome/genome.fa.fai').write_text('chr1\t6\t0\t6\t7\n')
        (source/'genome/assembly_report.txt').write_text('# Organism name: Example organism (fixture)\n')
        def gz(name,text):
            with gzip.open(source/name,'wt') as f:f.write(text)
        gz('genome/genome.fa.gz','>chr1\nacg\ntnn\n')
        gz('for_participants/whitelist.train.bed.gz','chr1\t1\t3\nchr1\t4\t6\n')
        gz('for_participants/template.train.bed.gz',''.join(f'chr1\t{a}\t{b}\t.\t0\t{s}\n' for s in ('+','-') for a,b in ((1,3),(4,6))))
        gz('for_participants/csRNA-train.txt.gz','0\n5\n0\n1\n7\n0\n0\n3\n')
        self.meta=build(self.raw,self.fixed,'jaNemVect1.1','train');(self.fixed/'catalog.json').write_text(json.dumps([self.meta]))
        self.parquet=self.root/'parquet';self.parquet.mkdir();self.record=export_one(self.fixed,self.parquet,self.meta,2)
        (self.parquet/'manifest.json').write_text(json.dumps({'files':[self.record]}))
    def tearDown(self):self.tmp.cleanup()
    def test_raw_to_per_base(self):
        rows=pq.read_table(self.parquet/self.record['file']).to_pylist()
        self.assertEqual([r['position'] for r in rows],list(range(6)))
        self.assertEqual(''.join(r['base'] for r in rows),'ACGTNN')
        self.assertEqual([r['csrna_plus'] for r in rows],[None,0,5,None,0,1])
        self.assertEqual([r['csrna_minus'] for r in rows],[None,7,0,None,0,3])
        self.assertEqual(self.record['eligible_bases'],4)
        self.assertTrue(self.record['all_row_groups_verified'])
    def test_vector_formats(self):
        path=self.root/'vector.gz'
        for text in ('0\n005\n0\n12', '0\r\n5\r\n0\r\n12\r\n', '0.0\n5.0\n0.0\n12.0\n'):
            with gzip.open(path,'wt') as stream:stream.write(text)
            pieces=list(vector_nonzero(path))
            self.assertEqual([int(v) for indices,_,_ in pieces for v in indices],[1,3])
            self.assertEqual([int(v) for _,counts,_ in pieces for v in counts],[5,12])
            self.assertEqual(sum(n for _,_,n in pieces),4)
        for text in ('-1\n', 'nan\n', '0\n\n1\n', '9007199254740993\n'):
            with gzip.open(path,'wt') as stream:stream.write(text)
            with self.assertRaises(ValueError):list(vector_nonzero(path))
    def test_downloads_and_ranges(self):
        original=(app.DATA,app.DOWNLOADS);app.DATA=self.fixed;app.DOWNLOADS=self.parquet
        for function in (app.catalog,app.track,app.reference_index,app.peaks):function.cache_clear()
        class Request:
            def __init__(self,path,method='GET',headers=''):
                self.input=io.BytesIO(f'{method} {path} HTTP/1.0\r\n{headers}\r\n'.encode());self.output=bytearray()
            def makefile(self,*a):return self.input
            def sendall(self,data):self.output.extend(data)
        class Quiet(app.Handler):
            def log_message(self,*a):pass
        def get(path,method='GET',headers=''):
            request=Request(path,method,headers);Quiet(request,('127.0.0.1',0),None)
            return bytes(request.output).split(b'\r\n\r\n',1)
        try:
            head,body=get('/api/downloads');self.assertEqual(len(json.loads(body)),1)
            path='/downloads/'+self.record['file'];expected=(self.parquet/self.record['file']).read_bytes()
            head,body=get(path,headers='Range: bytes=0-3\r\n');self.assertIn(b'206',head);self.assertEqual(body,b'PAR1')
            head,body=get(path,'HEAD');self.assertEqual(body,b'');self.assertIn(str(len(expected)).encode(),head)
            head,body=get(path,headers='Range: bytes=-4\r\n');self.assertEqual(body,b'PAR1')
            self.assertIn(b'416',get(path,headers='Range: bytes=999999-\r\n')[0])
            self.assertIn(b'404',get('/downloads/../../secret')[0])
            self.assertIn(b'404',get('/downloads/unknown.parquet')[0])
        finally:
            app.DATA,app.DOWNLOADS=original
            for function in (app.catalog,app.track,app.reference_index,app.peaks):function.cache_clear()

if __name__=='__main__':unittest.main()
