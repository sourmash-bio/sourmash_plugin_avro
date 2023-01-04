"Read/write sketches with Apache Avro."

# TODO:
# * add support for other moltypes

import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal

from sourmash.save_load import (Base_SaveSignaturesToLocation,
                                _get_signatures_from_rust)

import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

sig_schema_def = """
{
    "name": "Signature",
    "type":"record",
    "fields":[
       { "name": "class", "type": "string"},
       { "name": "email", "type": "string"},
       { "name": "hash_function", "type": "string"},
       { "name": "filename", "type": "string"},
       { "name": "name", "type": "string"},
       { "name": "license", "type": "string"},
       { "name": "version", "type": "float" },

       { "name": "signatures",
         "type": {
            "type": "array",
            "items": {

           "name": "MinHash",
           "type": "record",
           "fields":[
             { "name": "num", "type": "int" },
             { "name": "ksize", "type": "int" },
             { "name": "seed", "type": "int" },
             { "name": "max_hash", "type": { "name": "ulong", "type": "fixed", "size": 8 } },
             { "name": "md5sum", "type": "string" },
             { "name":"mins",
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
                     "type":"int"
                   }
                }
             }
           ]
         }
       }
       }
 
       ]
}
"""

schema = avro.schema.parse(sig_schema_def)

###

#
# load_from plugin - supports loading .avrosig files.
#

def load_sketches(location, *args, **kwargs):
    if location and location.endswith('.avrosig'):
        with open(location, "rb") as fp:
            reader = DataFileReader(fp, DatumReader())

            siglist = []
            for signature in reader:
                assert round(signature['version'], 1) == 0.4, signature['version']
                sketches = signature['signatures']
                moltype = signature['hash_function']

                for minhash in sketches:
                    if moltype.upper() == 'DNA': # @CTB
                        pass
                    else:
                        raise Exception

                    hashes = [ int.from_bytes(h, 'big') for h in minhash['mins'] ]
                    abunds = minhash['abunds']
                    is_abund = True
                    if max(abunds) == 1:
                        is_abund = False
                        abunds = None

                    max_hash = int.from_bytes(minhash['max_hash'], 'big')

                    mh = sourmash.MinHash(n=minhash['num'],
                                          ksize=minhash['ksize'],
                                          track_abundance=is_abund,
                                          max_hash=max_hash)

                    if is_abund:
                        abunds = dict(( (k, v) for k, v in zip(hashes,
                                                               abunds) ))
                        mh.set_abundances(abunds)
                    else:
                        mh.add_many(hashes)

                    ss = sourmash.SourmashSignature(mh,
                                                    name=signature['name'],
                                                    filename=signature['filename'])
                    siglist.append(ss)

        return LinearIndex(siglist)

load_sketches.priority = 5


#
# save_to plugin - supports output to .avrosig.
#

class SaveSignatures_AvroFile(Base_SaveSignaturesToLocation):
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

                hash_d = mh.hashes
                hashes = [ h.to_bytes(8, 'big') for h in hash_d ]
                abunds = [ hash_d[h] for h in hash_d ]

                minhash_d = dict(ksize=mh.ksize,
                                 num=mh.num,
                                 max_hash=mh._max_hash.to_bytes(8, 'big'),
                                 mins=hashes,
                                 abunds=abunds,
                                 seed=mh.seed,
                                 md5sum="")

                sig_d = dict(email='',
                             license='CC0',
                             hash_function="dna",
                             name=ss.name,
                             filename=ss.filename,
                             signatures=[minhash_d],
                             version=0.4)
                sig_d['class'] = 'sourmash_signature'

                writer.append(sig_d)
            
        writer.close()

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)
