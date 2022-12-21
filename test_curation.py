import pandas as pd
from Curation import Curation
import sys
import time

query = sys.argv[1]
strain = sys.argv[2]

boundaryFile = "Flu_profile_boundaries_Main.txt"
lookupTable = "Flu_profile_lookupTable_Main.txt"
profile_dir = "testing"

start = time.time()
result = Curation(query, strain, boundaryFile, lookupTable, profile_dir).curation_table()
end = time.time()

print("Time for one sequence:", round(end - start, 3), "seconds")
print('\n')
print(result)