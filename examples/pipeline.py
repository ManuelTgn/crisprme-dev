import crisprme2._crisprme2_native as native
import numpy as np

def myfunc(batch: native.PyAlignmentBatch):
    print(f"transform: {np.asarray(batch["offset"])}")

pipeline = native.create_pipeline(transform=myfunc)
pipeline.submit()
pipeline.receive()