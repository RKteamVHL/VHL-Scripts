import pandas as pd
from ..annotations.Annotation import AnnotationType


def keep_only_cases(df):
	out_df = df[df['type'] == AnnotationType.CASE]
	return out_df
