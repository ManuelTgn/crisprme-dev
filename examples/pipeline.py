import crisprme2._crisprme2_native as native
import numpy as np
import time

def myfunc(batch: native.PyAlignmentBatch):
    print(f"transform: {np.asarray(batch["offset"])}")

#class Scorer:
#    def __init__():
#        pass
#    def __call__(batch):
#        pass

#class Annotator:
#    def __init__():
#        pass
#    def __call__(batch):
#        pass

#scorer = Scorer(.....)
#pipeline = native.create_pipeline(transform=scorer, scores=[...], annotators=[...])
# add target batcher source
#pipeline.submit(batcher)
native.initialize_engine_logger()
pipeline = native.create_pipeline(transform=myfunc)
pipeline.submit_example() # make error visible
pipeline.close() # make error visible

complete = False
while not complete:
    complete, result = pipeline.receive() # make error visible
    if not complete:
        print(np.asarray(result["offset"]))
