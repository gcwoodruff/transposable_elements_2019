#this counts up the repeated items in a list

import numpy as np
import pandas as pd
import sys

inop_a = pd.read_csv(sys.argv[1])

inop_b = inop_a.apply(pd.value_counts)

inop_b.to_csv(sys.argv[2], sep='\t')
