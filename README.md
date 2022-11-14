
## Using modified mendel to calculate penetrance

### Program

1. The original version of the program was available at [30d8332](https://github.com/ictr/mendel_penetrance/tree/30d8332d8299e6ce4332c9bdf195bb8de82621a4/mendel)

2. We made a [single change](https://github.com/ictr/mendel_penetrance/commit/7491e344b8e27eb79a29a46c3fd9a5381818d235) so that it allows our input data.

3. We compile the code with

```
gfortran -O3 mendela.f single.f -o single
```

### Data

1. The test data is availabe in data directory.

2. We analyze the data using command

```
../mendel/single
```

and choose `NO`, `21`, and `no`

The output is in `single_out.dat`.


### References

* [AJHG 2003, Average risks of breast and ovarian cancer associated with BRCA1 or BRCA2 mutations](https://pubmed.ncbi.nlm.nih.gov/12677558/)
* [Br J Cancer, 2008, The BOADICEA model of genetic susceptibility to breast and ovarian cancers](https://pubmed.ncbi.nlm.nih.gov/18349832/)
* [NEJM, 2014, Breast-cancer risk in families with mutations in PALB2](https://pubmed.ncbi.nlm.nih.gov/25099575/)
* [J Clin Onco, 2020, Cancer Risks Associated With Germline PALB2 Pathogenic Variants](https://pubmed.ncbi.nlm.nih.gov/31841383/)
