import crisprme2._crisprme2_native as native
import numpy as np
import time

def myfunc(batch: native.PyAlignmentBatch):
    print(f"transform: {np.asarray(batch["offset"])}")

pipeline = native.create_pipeline(transform=myfunc)
pipeline.submit() # add target batcher source
pipeline.close()

complete = False
while not complete:
    complete, result = pipeline.receive()
    if not complete:
        print(np.asarray(result["offset"]))