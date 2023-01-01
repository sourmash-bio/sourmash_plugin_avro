"Read/write sketches with Apache Avro."

# TODO:
# * add support for num
# * add support for other moltypes
# * add support for abundances

import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal

# CTB: remove _ prefix?
from sourmash.sourmash_args import _BaseSaveSignaturesToLocation
from sourmash.sourmash_args import _get_signatures_from_rust

import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

minhash_schema_def = """{
    "namespace": "bio.sourmash.avro.schema",
    "name": "MinHash",
    "type":"record",
    "fields":[
        { "name": "num", "type": "int" },
        { "name": "scaled", "type": "int" },
        { "name": "ksize", "type": "int" },
        { "name": "molecule", "type": {
          "type": "enum",
          "name": "moltype",
          "symbols": ["DNA", "protein", "hp", "dayhoff"] }
        },
        { "name":"hashes",
          "type": {
             "type": "array",  
              "items":{
                  "name":"hash",
                   "type":"fixed",
                   "size": 8
              }
           }
        },
        { "name":"abunds",
          "type": {
             "type": "array",  
              "items":{
                  "name":"abund",
                   "type": "int"
              }
           }
        }
    ]
    }"""

sig_schema_def = """{
    "namespace": "bio.sourmash.avro.schema",
    "name": "SourmashSignature",
    "type":"record",
    "fields":[
       { "name": "filename", "type": "string"},
       { "name": "name", "type": "string"}
     ]
    }
"""
sig_schema_def = """{
    "namespace": "bio.sourmash.avro.schema",
    "name": "SourmashSignature",
    "type":"record",
    "fields":[
       { "name": "filename", "type": "string"},
       { "name": "name", "type": "string"},
       { "name": "minhash", "type": "bio.sourmash.avro.schema.MinHash" }
     ]
    }
"""

#minhash_schema = avro.schema.parse(minhash_schema_def)

combined = minhash_schema_def + ",\n" + sig_schema_def
print('XXX', combined[800:870])
schema = avro.schema.parse(minhash_schema_def + ",\n" + sig_schema_def)
#schema = avro.schema.parse(minhash_schema_def, sig_schema_def)
#sketch_schema = avro.schema.parse(minhash_schema_def)

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
                    
                mh = sourmash.MinHash(n=sketch['num'],
                                      ksize=sketch['ksize'],
                                      scaled=sketch['scaled'])

                hashes = [ int.from_bytes(h, 'big') for h in sketch['hashes'] ]
                abunds = [ 0 for h in hashes ]
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
                                DatumWriter(), schema)

        for ss0 in self.keep:
            for ss in _get_signatures_from_rust([ss0]):
                mh = ss.minhash
                hashes = [ h.to_bytes(8, 'big') for h in mh.hashes ]
                abunds = [ 0 for h in hashes ]

                minhash_d = dict(ksize=mh.ksize,
                                           num=mh.num,
                                           scaled=mh.scaled,
                                           hashes=hashes,
                                           abunds=abunds,
                                           molecule=mh.moltype)
                writer.append(dict(minhash=minhash_d,
                                   name=ss.name,
                                   filename=ss.filename))
            
        writer.close()

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)
