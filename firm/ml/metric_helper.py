from __future__ import absolute_import, division, print_function

import numpy as np

from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score, accuracy_score

def get_formatted_scores(y_true, y_pred, y_pred_prob):
    return format_scores(calc_scores(y_true, y_pred, y_pred_prob))

def calc_scores(y_true, y_pred, y_pred_prob):
    return np.array([metric_func(y_true, y_pred_prob) if is_prob_metric else metric_func(y_true, y_pred) for \
            metric_func, is_prob_metric in _METRIC_FUNCS_AND_MODES], dtype = np.float64)

def format_scores(scores):
    return ', '.join(['%s = %f' % metric_tuple for metric_tuple in zip(_METRIC_NAMES, scores)])

def specificity_score(y_true, y_pred, *args, **kwargs):
    return recall_score(y_true == 0, y_pred == 0, *args, **kwargs)

_METRICS = [
    ('AUC', (roc_auc_score, True)),
    ('F1', (f1_score, False)),
    ('Precision', (precision_score, False)),
    ('Recall', (recall_score, False)),
    ('Specificity', (specificity_score, False)),
    ('Accuracy', (accuracy_score, False)),
]

_METRIC_NAMES, _METRIC_FUNCS_AND_MODES = zip(*_METRICS)

