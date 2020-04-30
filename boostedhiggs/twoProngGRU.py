import tensorflow as tf
from tensorflow import keras
from coffea.nanoaod import NanoEvents
from coffea.nanoaod.methods import Candidate
import numpy as np
import awkward
import os

# let tf us only a single thread
os.environ["OMP_NUM_THREADS"] = "1"
tf.config.threading.set_inter_op_parallelism_threads(1)
tf.config.threading.set_intra_op_parallelism_threads(1)

model = keras.models.load_model("../boostedhiggs/data/weights_gru.h5")

def run_model(jets):
    def pad_and_flatten(val): 
        print (val)
        return val.pad(20, clip=True).fillna(0.).flatten().reshape(-1, 20)
    inputs = np.stack([

        # normalize pt, eta, phi of cands to jet axis
        pad_and_flatten(jets.PFCands.pt / jets.pt),
        pad_and_flatten(jets.PFCands.eta - jets.eta),
        pad_and_flatten(jets.PFCands.delta_phi(jets)),
   
        # vertex params are already regular
        pad_and_flatten(jets.PFCands.dz),
        pad_and_flatten(jets.PFCands.d0),
    ])
    inputs = np.moveaxis(inputs, 0, 2)
    outputs = model.predict(inputs)

    # 2-pronginess is second column of output
    twoprong = outputs[:, 1]
    return twoprong
