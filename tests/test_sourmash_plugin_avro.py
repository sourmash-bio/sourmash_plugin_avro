"""
Tests for the sourmash avro read/write plugin.

TODO:
* test num
"""
import os
import pytest

import avro
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_sourmash_save_avro(runtmp):
    # make sure .avrosig => avro file.
    sig47 = utils.get_test_data('47.fa.sig')
    avro_out = runtmp.output('47.avrosig')

    runtmp.sourmash('sig', 'cat', sig47, '-o', avro_out)

    with open(avro_out, 'rb') as fp:
        reader = DataFileReader(fp, DatumReader())

        sigs = list(reader)
        assert len(sigs) == 1

        ss = sigs[0]
        assert len(ss['signatures']) == 1

        minhash = ss['signatures'][0]
        assert minhash['ksize'] == 31 # just a lite check ;)


def test_sourmash_load_avro(runtmp):
    # can we load w/sourmash? check object equality.
    sig47 = utils.get_test_data('47.fa.sig')
    avro_out = utils.get_test_data('47.avrosig')

    # ok, load using sourmash and compare.
    json_sig = sourmash.load_file_as_index(sig47)
    json_sig = list(json_sig.signatures())[0]

    avro_sig = sourmash.load_file_as_index(avro_out)
    avro_sig = list(avro_sig.signatures())[0]

    assert json_sig == avro_sig


def test_sourmash_load_avro_describe(runtmp):
    # can we load w/sourmash? check output of 'describe'
    sig47_avro = utils.get_test_data('47.avrosig')

    runtmp.sourmash('sig', 'describe', sig47_avro)

    out = runtmp.last_result.out
    err = runtmp.last_result.out

    print(out, err)

    expected = """signature: NC_009665.1 Shewanella baltica OS185, complete genome
source file: 47.fa
md5: 09a08691ce52952152f0e866a59f6261
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 5177
sum hashes: 5177
signature license: CC0""".splitlines()

    for line in expected:
        assert line in out


def test_sourmash_save_avro_abund(runtmp):
    # check abundance stuff
    sig47 = utils.get_test_data('abund/47.fa.sig')
    avro_out = runtmp.output('47.abund.avrosig')

    runtmp.sourmash('sig', 'cat', sig47, '-o', avro_out)

    # legit avro?
    with open(avro_out, 'rb') as fp:
        reader = DataFileReader(fp, DatumReader())

        sigs = list(reader)

    # ok, load using sourmash and compare.
    json_sig = sourmash.load_file_as_index(sig47)
    json_sig = list(json_sig.signatures())[0]

    avro_sig = sourmash.load_file_as_index(avro_out)
    avro_sig = list(avro_sig.signatures())[0]

    assert json_sig == avro_sig


def test_sourmash_load_save_prot(runtmp):
    # check abundance stuff
    sigfile = utils.get_test_data('prot.sig.gz')
    avro_out = runtmp.output('out.avrosig')

    runtmp.sourmash('sig', 'cat', sigfile, '-o', avro_out)

    # legit avro?
    with open(avro_out, 'rb') as fp:
        reader = DataFileReader(fp, DatumReader())

        sigs = list(reader)

    # ok, load using sourmash and compare.
    json_sig = sourmash.load_file_as_index(sigfile)
    json_sig = list(json_sig.signatures())[0]

    avro_sig = sourmash.load_file_as_index(avro_out)
    avro_sig = list(avro_sig.signatures())[0]

    assert json_sig == avro_sig


def test_sourmash_load_save_dayhoff(runtmp):
    # check abundance stuff
    sigfile = utils.get_test_data('dayhoff.sig.gz')
    avro_out = runtmp.output('out.avrosig')

    runtmp.sourmash('sig', 'cat', sigfile, '-o', avro_out)

    # legit avro?
    with open(avro_out, 'rb') as fp:
        reader = DataFileReader(fp, DatumReader())

        sigs = list(reader)

    # ok, load using sourmash and compare.
    json_sig = sourmash.load_file_as_index(sigfile)
    json_sig = list(json_sig.signatures())[0]

    avro_sig = sourmash.load_file_as_index(avro_out)
    avro_sig = list(avro_sig.signatures())[0]

    assert json_sig == avro_sig


def test_sourmash_load_save_hp(runtmp):
    # check abundance stuff
    sigfile = utils.get_test_data('hp.sig.gz')
    avro_out = runtmp.output('out.avrosig')

    runtmp.sourmash('sig', 'cat', sigfile, '-o', avro_out)

    # legit avro?
    with open(avro_out, 'rb') as fp:
        reader = DataFileReader(fp, DatumReader())

        sigs = list(reader)

    # ok, load using sourmash and compare.
    json_sig = sourmash.load_file_as_index(sigfile)
    json_sig = list(json_sig.signatures())[0]

    avro_sig = sourmash.load_file_as_index(avro_out)
    avro_sig = list(avro_sig.signatures())[0]

    assert json_sig == avro_sig
