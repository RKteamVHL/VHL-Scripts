import os

import graphviz
from sklearn import tree
from sklearn.tree import DecisionTreeRegressor

from . import constants
from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS


def dataframe_decisiontree(df, analysis_type, missense_only=False):
    if missense_only:
        df = df.dropna(subset=['generalized_mutant_type.missense_variant'])
    fig_path = os.path.join(constants.FIGURE_DIR, analysis_type)
    df = df.dropna(subset=['codon_start'])
    df = df.fillna(False)
    df[COMPUTED_COLUMNS["generalized_mutant_type"]] = df[COMPUTED_COLUMNS["generalized_mutant_type"]].astype(bool)
    cols = ['codon_start', *COMPUTED_COLUMNS["generalized_mutant_type"]]
    X = df[cols]
    y = df[COMPUTED_COLUMNS["generalized_phenotype"]]
    for i in range(1, 5):

        regr = DecisionTreeRegressor(max_depth=i).fit(X, y)

        dot_data = tree.export_graphviz(regr, out_file=None, feature_names=cols, class_names=True, filled=True, rounded=True)
        graph = graphviz.Source(dot_data)
        graph.render(os.path.join(fig_path, f"depth_{i}"), cleanup=True)
    regr = DecisionTreeRegressor().fit(X, y)
    dot_data = tree.export_graphviz(regr, out_file=None, feature_names=cols, class_names=True, filled=True,
                                    rounded=True)
    graph = graphviz.Source(dot_data)
    graph.render(os.path.join(fig_path, f"depth_none"), cleanup=True)
    pass