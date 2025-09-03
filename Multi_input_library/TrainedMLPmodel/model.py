#    This file was created by
#    MATLAB Deep Learning Toolbox Converter for TensorFlow Models.
#    03-Sep-2025 11:27:13

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

def create_model():
    input = keras.Input(shape=(41,))
    fc_1 = layers.Dense(160, name="fc_1_")(input)
    layer_1 = layers.Activation('tanh')(fc_1)
    fc_2 = layers.Dense(80, name="fc_2_")(layer_1)
    layer_2 = layers.Activation('tanh')(fc_2)
    fc_3 = layers.Dense(40, name="fc_3_")(layer_2)
    layer_3 = layers.Activation('tanh')(fc_3)
    fc_4 = layers.Dense(20, name="fc_4_")(layer_3)
    layer_4 = layers.Activation('tanh')(fc_4)
    fc_5 = layers.Dense(10, name="fc_5_")(layer_4)
    layer_5 = layers.Activation('tanh')(fc_5)
    fc_6 = layers.Dense(4, name="fc_6_")(layer_5)
    relu = layers.ReLU()(fc_6)

    model = keras.Model(inputs=[input], outputs=[relu])
    return model
