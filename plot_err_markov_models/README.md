## Age Buckets

Age buckets are in 5 year increments from 0 to 30 (6 age buckets).

## Data Format

To open:

```python
import numpy as np
d = np.load('CPG_DICT_TO_PLOT.npy', allow_pickle='TRUE').item()
```

```python
d
>>> {'cg16984944': {0: {(-1.04, -0.96): 0.15223609624255596,
   (-0.96, -0.8799999999999999): 0.16293458764934712,
   (-0.8799999999999999, -0.7999999999999998): 0.183840791959947,
   (-0.31999999999999984, -0.23999999999999977): 0.5009885241481501},
  1: {(0.4800000000000004, 0.5600000000000005): 0.39037538238864833,
   (1.4400000000000004, 1.5200000000000005): 0.13465720951355067,
   (0.3200000000000003, 0.40000000000000036): -0.5432918133756863},
  2: {(-0.3999999999999999, -0.31999999999999984): 0.127078202461335},
  3: {(-0.08000000000000007, 0.0): 0.5522798595303037,
   (0.2400000000000002, 0.3200000000000003): -1.0},
  4: {(-0.6400000000000001, -0.56): 0.10707113430608196,
   (-0.31999999999999984, -0.23999999999999977): 0.23564429573263196,
   (0.0, 0.08000000000000007): -0.7821861858431841,
   (0.2400000000000002, 0.3200000000000003): -0.21781381415681583}},
 'cg24046474': {0: {(0.8799999999999999, 0.96): 0.31095712565028866,
   (1.6799999999999997, 1.7599999999999998): 0.1687047359487461,
   ...
```
And,
```python
d.keys()
>>> ['cg16984944',
 'cg24046474',
 'cg15361750',
 'cg17274064',
 'cg07979752',
 'cg10637955',
 'cg26079320',
 'cg13150977',
 'cg26149738',
 'cg17142470']
```

## About the Dictionary

Here, the dictionary is indexed in the following order:
1. CpG site
2. Transition (e.g. 0 indicates a transition from the (0-5) age group to the (5-10) age group)
3. Methylation State

Value: Transition probability error (healthy transition probability - error probability)

# What needs to be done
All that needs to be done is for each CpG site, create a markov model visualization of the data we have here to go in the final paper.

So, that's 10 hidden markov models.

It would be great if you could color the negative transition weights a different color that positive ones.
