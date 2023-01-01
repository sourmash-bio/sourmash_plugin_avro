"Read/write sketches with Apache Avro."
import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal

# CTB: remove _ prefix?
from sourmash.sourmash_args import _BaseSaveSignaturesToLocation
from sourmash.sourmash_args import _get_signatures_from_rust

import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

sketch_schema_def = """
{
    "name": "Sketch",
    "type":"record",
    "fields":[
        { "name": "scaled", "type": "int" },
        { "name": "ksize", "type": "int" },
        { "name": "filename", "type": "string"},
        { "name": "name", "type": "string"},
        { "name": "molecule", "type": {
            "type": "enum",
            "name": "moltype",
            "symbols": ["DNA", "protein", "hp", "dayhoff"] }},
       { "name":"hashes",
            "type": {
                "type": "array",  
                "items":{
                    "name":"hash",
                    "type":"fixed",
                    "size": 8
                }
            }
        }
    ] 
}
"""
sketch_schema = avro.schema.parse(sketch_schema_def)


###

def load_sketches(location, *args, **kwargs):
    if location and location.endswith('.avrosig'):
        with open(location, "rb") as fp:
            reader = DataFileReader(fp, DatumReader())

            siglist = []
            for sketch in reader:
                moltype = sketch['molecule']
                if moltype == 'DNA':
                    pass
                else:
                    raise Exception
                    
                mh = sourmash.MinHash(n=0,
                                      ksize=sketch['ksize'],
                                      scaled=sketch['scaled'])

                hashes = [ int.from_bytes(h, 'big') for h in sketch['hashes'] ]
                mh.add_many(hashes)

                ss = sourmash.SourmashSignature(mh,
                                                name=sketch['name'],
                                                filename=sketch['filename'])
                siglist.append(ss)

        return LinearIndex(siglist)

load_sketches.priority = 5


class SaveSignatures_AvroFile(_BaseSaveSignaturesToLocation):
    "Save signatures to an Apache Avro file." 
    def __init__(self, location):
        super().__init__(location)
        self.keep = []          # CTB: could switch to open file.

    @classmethod
    def matches(self, location):
        # match anything that is not None or ""
        if location:
            return location.endswith('.avrosig')

    def __repr__(self):
        return f"SaveSignatures_AvroFile('{self.location}')"

    def open(self):
        pass

    def close(self):
        writer = DataFileWriter(open(self.location, "wb"),
                                DatumWriter(), sketch_schema)

        for ss0 in self.keep:
            for ss in _get_signatures_from_rust([ss0]):
                mh = ss.minhash
                ksize = mh.ksize
                scaled = mh.scaled
                hashes = [ h.to_bytes(8, 'big') for h in mh.hashes ]
                name = ss.name
                filename = ss.filename
                molecule = mh.moltype

                writer.append(dict(ksize=ksize,
                                   scaled=scaled,
                                   hashes=hashes,
                                   name=name,
                                   filename=filename,
                                   molecule=molecule))
            
        writer.close()

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)
