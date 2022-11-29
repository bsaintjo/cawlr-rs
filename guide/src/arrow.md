# Reading and Writing Arrow files from cawlr

## Scored Reads files

### Format

For smaller datasets and enough ram, you can use Python dictionaries to convert them into Arrow files

```python
{ 'scored': [ { 'metadata': { 'chrom': 'chrXIII',
                              'length': 178,
                              'name': '20d1aac0-29de-43ae-a0ef-aa8a6766eb70',
                              'seq': '',
                              'start': 182504,
                              'strand': True},
                'scores': [ { 'kmer': 'TATTCA',
                              'pos': 182509,
                              'score': 0.898015077423625,
                              'signal_score': 0.898015077423625,
                              'skip_score': 0.0,
                              'skipped': False},
                            { 'kmer': 'ATCCTA',
                              'pos': 182676,
                              'score': 1.0,
                              'signal_score': 1.0,
                              'skip_score': 0.0,
                              'skipped': False}]}]}
```