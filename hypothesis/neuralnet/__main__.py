from ..config import OUTPUT_DIR
from ..variant_functions import VHL_PHENOTYPES
import os
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers

RNG = np.random.default_rng()

VHL_AA_LENGTH = 213

AA_LIST = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
           'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'TER']

def df_to_dataset(dataframe, shuffle=True, batch_size=32):
    df = dataframe.copy()
    labels = df.pop('target')
    df = {key: value[:, tf.newaxis] for key, value in dataframe.items()}
    ds = tf.data.Dataset.from_tensor_slices((dict(df), labels))
    if shuffle:
        ds = ds.shuffle(buffer_size=len(dataframe))
    ds = ds.batch(batch_size)
    ds = ds.prefetch(batch_size)
    return ds


if __name__ == '__main__':

    missense_df = pd.read_csv(os.path.join(OUTPUT_DIR, 'missense_variants.csv'))
    missense_df = missense_df.drop(['Unnamed: 0'], axis=1)

    missense_df.loc[:, 'pos'] = missense_df['pos'].astype(np.float) / VHL_AA_LENGTH

    train, val, test = np.split(missense_df.sample(frac=1), [int(0.8 * len(missense_df)), int(0.9 * len(missense_df))])

    print(len(train), 'training examples')
    print(len(val), 'validation examples')
    print(len(test), 'test examples')


    # # train: 50%, validate: 25%, test: 25%
    #
    # train_df = missense_df.sample(frac=0.5, random_state=RNG)
    #
    # valid_df = missense_df.drop(train_df.index).sample(frac=0.5, random_state=RNG)
    #
    # test_df = missense_df.drop(train_df.index.union(valid_df.index))
    #
    # # sanity check of data splitting- should all be 0
    # print(len(train_df.index.intersection(valid_df.index)))
    # print(len(train_df.index.intersection(test_df.index)))
    # print(len(valid_df.index.intersection(test_df.index)))


    # preprocessing for codons, binning
    resolution_in_codons = 3
    boundaries = list(np.arange(0, VHL_AA_LENGTH+resolution_in_codons+1, resolution_in_codons))

    discretization_layer = tf.keras.layers.Discretization(bin_boundaries=boundaries)
    position_one_hot_layer = tf.keras.layers.CategoryEncoding(
        num_tokens=len(boundaries)+1, output_mode='one_hot')
    position_one_hot_layer(discretization_layer([-5, 0, 3, 100, 200, 213, 500]))

    # preprocessing for amino acid names
    string_lookup_layer = tf.keras.layers.StringLookup(
        vocabulary=AA_LIST,
        num_oov_indices=0,
        output_mode='one_hot')
    string_lookup_layer(['TER', 'MET', 'ALA'])

    # preprocessing for from_aa-to_aa feature cross
    cross_layer = tf.keras.layers.experimental.preprocessing.HashedCrossing(
    num_bins=len(AA_LIST)**2, output_mode='one_hot')


    batch_size = 256
    train_ds = df_to_dataset(train, batch_size=batch_size)
    val_ds = df_to_dataset(val, shuffle=False, batch_size=batch_size)
    test_ds = df_to_dataset(test, shuffle=False, batch_size=batch_size)

    all_inputs = []
    encoded_features = []

    # position input
    pos_col = tf.keras.Input(shape=(1,), name='Position', dtype='int64')
    encoding_layer = position_one_hot_layer
    encoded_pos_col = encoding_layer(pos_col)
    all_inputs.append(pos_col)
    encoded_features.append(encoded_pos_col)
