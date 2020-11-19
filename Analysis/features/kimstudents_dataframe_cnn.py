import os
from . import constants

def dataframe_decisiontree(df, analysis_type):
    fig_path = os.path.join(constants.FIGURE_DIR, analysis_type)
    df = df.dropna(subset=['codon_start']).fillna(0)