# Pickle file format

The outputs from `cawlr train` and `cawlr rank` use the Python pickle format. You can use the [pickle library](https://docs.python.org/3/library/pickle.html).

## Loading the file

```python
import pickle

# Can also be output from cawlr train
with open("ranks.pickle", 'rb') as rank_file:
    ranks = pickle.load(rank_file)
```

## `cawlr rank`

The `ranks` object acts as a Python dictionary, mapping a kmer to the Kulback-Liebler divergence between the models from the positive and negative controls.

```python
>>> ranks
{'TACTAC': 0.36897146906453593, 'GCTGAC': 1.1134219577325406, 'TCTGCG': 0.7802014312171258, ..
# cutoff for brevity
}

```

## `cawlr train`

Similar to above, the loaded file acts like a Python dictionary mapping kmers to

```python
>>> with open("/scratch/bsaintjo/220525_unique_cawlr/unique.500.model", 'rb') as model_file:
...     model = pickle.load(model_file)
```
