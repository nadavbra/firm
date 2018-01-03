from __future__ import absolute_import, division, print_function

import numpy as np

from sklearn.cross_validation import KFold
from sklearn.metrics import confusion_matrix

from .util import log
from .classification import predict_prob, find_prob_threshold
from .metric_helper import get_formatted_scores, calc_scores, format_scores

def cross_validate(X, y, model, n_folds = 3, feature_selection = None, feature_names = None, report_removed_features = False, seed = None):
    
    log('Starting %d folds of cross validation...\n' % n_folds)
    kf = KFold(len(y), shuffle = True, random_state = seed, n_folds = n_folds)
    all_test_scores = []

    for train_index, test_index in kf:
        
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        log('KFold train_size = %d, test_size = %d' % (len(y_train), len(y_test)))
        
        test_scores = estimate_model(X_train, y_train, X_test, y_test, model, feature_selection = feature_selection, \
                feature_names = feature_names, report_removed_features = report_removed_features)
        all_test_scores.append(test_scores)
    
        log('*' * 50)
        
    avg_test_scores = np.mean(np.array(all_test_scores), axis = 0)
    log('Average test scores: ' + format_scores(avg_test_scores))
    log('Done.')
        
def estimate_model(X_train, y_train, X_test, y_test, model, feature_selection = None, feature_names = None, report_removed_features = False):

    if feature_selection is not None:
        X_train, X_test = select_features(feature_selection, feature_names, X_train, y_train, X_test, report_removed_features = False)
    
    prob_threshold = fit_model_and_report_training_performance(X_train, y_train, model)
    return report_and_collect_test_scores(X_test, y_test, model, prob_threshold)
    
def select_features(feature_selection, feature_names, X_train, y_train, X_test, report_removed_features = False):
    
    feature_selection.fit(X_train, y_train)
    removed_features = feature_names[~feature_selection.get_support()]
    
    if report_removed_features:
        log('Filtered out the following %d features: %s' % (len(removed_features), ', '.join(removed_features)))
    else:
        log('Filtered out %d features.' % len(removed_features))
    
    X_train = feature_selection.transform(X_train)
    X_test = feature_selection.transform(X_test)
    
    log('Selected features: %d' % X_train.shape[1])
    
    return X_train, X_test
    
def fit_model_and_report_training_performance(X_train, y_train, model):
    
    model.fit(X_train, y_train)
    
    y_pred_prob_train = predict_prob(model, X_train)
    prob_threshold = find_prob_threshold(y_train, y_pred_prob_train)
    log('Probability threshold: %f' % prob_threshold)
    
    y_pred_train = (y_pred_prob_train >= prob_threshold)
    
    log('Training confusion matrix:')
    log(confusion_matrix(y_train, y_pred_train))
    
    log('Training scores: ' + get_formatted_scores(y_train, y_pred_train, y_pred_prob_train))
    
    return prob_threshold
    
def report_and_collect_test_scores(X_test, y_test, trained_model, prob_threshold = 0.5):

    y_pred_prob_test = predict_prob(trained_model, X_test)
    y_pred_test = (y_pred_prob_test >= prob_threshold)
    
    log('Test confusion matrix:')
    log(confusion_matrix(y_test, y_pred_test))
    
    test_scores = calc_scores(y_test, y_pred_test, y_pred_prob_test)
    log('Test scores: ' + format_scores(test_scores))
    return test_scores
    