import pandas as pd
from Curation import Curation
import sys
import time

query = sys.argv[1]

start = time.time()
result = Curation(query).curation_table()
end = time.time()

print("Time for one sequence:", round(end - start, 3), "seconds")
print('\n')
print(result)