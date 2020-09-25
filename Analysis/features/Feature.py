import logging

import numpy as np


class Feature:
    def __init__(self, *args, name="Feature", column=None, **kwargs):
        self.name = name
        self.column = column
        self.rows = {}
        self.logger = logging.getLogger(self.name)

    @staticmethod
    def validate_row(row):
        if row is not None:
            return True

    def update_by_column(self, row, column):
        category = row[column].casefold()
        self.update(row, category)

    def add_category(self, category):
        self.rows[category] = self.rows.get(category, [])

    def update(self, row, category):
        if self.validate_row(row):
            self.add_category(category)
            self.rows[category].append(row)


class NominalFeature(Feature):
    # sorts row by category counts
    def sort(self):
        keys = np.array(list(self.rows.keys()))
        counts = [len(v) for v in self.rows.values()]
        c_order = np.argsort(counts)
        k_order = keys[c_order]
        old_rows = self.rows
        self.rows = {k: old_rows[k] for k in k_order}

    def to_dict(self):
        return {k: len(self.rows[k]) for k in self.rows.keys()}


class OrdinalFeature(NominalFeature):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)


class IntervalFeature(Feature):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.histogram = OrdinalFeature(name=f'{self.name}_hist')

    def add_value(self, value):
        # TODO: some error checking to account for value rounding issues
        self.rows[value] = self.rows.get(value, [])

    def update(self, row, value):
        if self.validate_row(row):
            self.add_value(value)
            self.rows[value].append(row)

    def make_hist(self, bins=10):

        values = np.array([k for k in self.rows.keys()])
        bins_out = np.histogram_bin_edges(values, bins=bins)
        values_bins = np.digitize(values, bins_out)

        for i in range(len(bins_out) - 1):
            category = f'{int(bins_out[i])}-{int(bins_out[i + 1])}'
            self.histogram.add_category(category)

        for i in range(len(values)):
            v_bin = values_bins[i]
            if v_bin == len(bins_out):
                v_bin -= 1
            bin_lower = int(bins_out[v_bin - 1])
            bin_upper = int(bins_out[v_bin])

            for row in self.rows[values[i]]:
                self.histogram.update(row, f'{bin_lower}-{bin_upper}')

    def sort(self):
        keys = np.array(list(self.rows.keys()))
        counts = [float(k) for k in self.rows.keys()]
        c_order = np.argsort(counts)
        k_order = keys[c_order]
        old_rows = self.rows
        self.rows = {k: old_rows[k] for k in k_order}

    def to_dict(self, bins=None):
        if bins is not None:
            self.make_hist(bins=bins)
            out_dict = self.histogram.to_dict()

        else:
            self.sort()
            out_dict = {k: len(self.rows[k]) for k in self.rows.keys()}

        return out_dict


class RatioFeature(IntervalFeature):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
