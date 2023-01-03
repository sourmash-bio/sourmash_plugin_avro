# sourmash signature/sketch storage in the Apache Avro format.

This is an **experimental** plugin to support Apache Avro storage
of sourmash sketches.

Links:
* [sourmash: sourmash-bio/sourmash/](https://github.com/sourmash-bio/sourmash/)
* [Apache Avro](https://avro.apache.org/)

Note: Depends on
[sourmash#2428](https://github.com/sourmash-bio/sourmash/pull/2428).

## Usage

```
pip install .
```

and then you can create Avro output files by using `-o filename.avrosig`.
These files can be loaded with sourmash without any further work.

