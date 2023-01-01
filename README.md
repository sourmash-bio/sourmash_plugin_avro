# Apache Avro plugin supporting sourmash signature/sketch storage

This is an **experimental** plugin to support Apache Avro storage
of sourmash sketches.

[sourmash: sourmash-bio/sourmash/](https://github.com/sourmash-bio/sourmash/)

[Avro](https://avro.apache.org/)

Note: Depends on
[sourmash#2428](https://github.com/sourmash-bio/sourmash/pull/2428).

## Usage

```
pip install .
```

and then you can create Avro output files by using `-o filename.avrosig`.
These files can be loaded with sourmash without any further work.
