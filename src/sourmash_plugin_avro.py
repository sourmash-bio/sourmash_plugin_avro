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

sig_schema_def = """
...
"""

# load from file for easier editing... at least for now!
schema = avro.schema.parse(open('signature.avsc', 'rb').read())

###

def load_sketches(location, *args, **kwargs):
    if location and location.endswith('.avrosig'):
        with open(location, "rb") as fp:
            reader = DataFileReader(fp, DatumReader())

            siglist = []
            for signature in reader:
                sketches = signature['signatures']

                for minhash in sketches:
                    moltype = minhash['molecule']
                    if moltype == 'DNA': # @CTB
                        pass
                    else:
                        raise Exception

                    hashes = [ int.from_bytes(h, 'big') for h in minhash['mins'] ]
                    abunds = minhash['abunds']
                    is_abund = True
                    if max(abunds) == 0:
                        is_abund = False
                        abunds = None

                    mh = sourmash.MinHash(n=minhash['num'],
                                          ksize=minhash['ksize'],
                                          track_abundance=is_abund,
                                          scaled=1000) # @CTB

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

                hash_d = mh.hashes
                hashes = [ h.to_bytes(8, 'big') for h in hash_d ]
                abunds = [ hash_d[h] for h in hash_d ]

                minhash_d = dict(ksize=mh.ksize,
                                 num=mh.num,
#           { "name": "max_hash", "type": "fixed", "size": 8 },
 #                                 max_hash=mh._max_hash.to_bytes(8, 'big'),
                                 mins=hashes,
                                 abunds=abunds,
                                 molecule=mh.moltype,
                                 seed=mh.seed,
                                 md5sum="")

                sig_d = dict(email='',
                             hash_function='0.murmur64',
                             license='CC0',
                             name=ss.name,
                             filename=ss.filename,
                             signatures=[minhash_d])
                sig_d['class'] = 'sourmash_signature'

                writer.append(sig_d)
            
        writer.close()

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)
