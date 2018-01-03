from __future__ import absolute_import, division, print_function

import pickle

import numpy as np

from sklearn.metrics import recall_score

from .util import log, is_iterable, split_to_size, create_random_mask, distirbute_with_callback
from .feature_engineering import invoke_feature_extractor_to_matrix
from .metric_helper import specificity_score

class Classifier(object):

    def __init__(self, feature_extractor_creator, feature_extractor_creator_args = [], feature_extractor_creator_kwargs = {}, \
            model = None, feature_selection = None, default_prob_threshold = None, prob_adjustment_table = None):

        self.feature_extractor_creator = feature_extractor_creator
        self.feature_extractor_creator_args = feature_extractor_creator_args
        self.feature_extractor_creator_kwargs = feature_extractor_creator_kwargs
        self.feature_extractor = feature_extractor_creator(*feature_extractor_creator_args, **feature_extractor_creator_kwargs)
        self.model = model
        self.feature_selection = feature_selection
        self.default_prob_threshold = default_prob_threshold
        self.prob_adjustment_table = prob_adjustment_table
    
    def predict_proba(self, sample_or_samples, chunk_size = 512, thread_pool = None):

        bulk_input = is_iterable(sample_or_samples)
        
        if bulk_input:
            samples = np.array(sample_or_samples)
        else:
            samples = np.array([sample_or_samples])
            
        sample_chunks = split_to_size(samples, chunk_size)
        y_pred = self._predict_proba_on_chunks(sample_chunks, thread_pool)
    
        if bulk_input:
            return y_pred
        else:
            return y_pred[0]
    
    def predict_adjusted_proba(self, sample_or_samples, **kwargs):
        return np.interp(self.predict_proba(sample_or_samples, **kwargs), *zip(*self.prob_adjustment_table))
        
    def predict(self, sample_or_samples, prob_threshold = None, **kwargs):
        
        if prob_threshold is None:
            prob_threshold = self.default_prob_threshold
            
        return self.predict_proba(sample_or_samples, **kwargs) >= prob_threshold
        
    def train_on_processed_data(self, X, y, learn_prob_adjustment = True, prob_adjustment_bins = 20):
    
        if self.default_prob_threshold is not None:
            log('Warning: %s was already set with a probability threshold; it will be overriden.' % self)
    
        if self.feature_selection is not None:
            self.feature_selection.fit(X, y)
            X = self.feature_selection.transform(X)
        
        if learn_prob_adjustment:
            self._learn_prob_adjustment_table(X, y, prob_adjustment_bins)
        
        self.model.fit(X, y)
        
        y_pred_prob = predict_prob(self.model, X)
        self.default_prob_threshold = find_prob_threshold(y, y_pred_prob)
        
    def train_on_raw_data(self, samples, labels):
        X = self.feature_extractor.create_features_as_matrix(samples)
        self.train_on_processed_data(X, labels)
        
    def dump(self, file_handler):
        state = (self.model, self.feature_selection, self.default_prob_threshold, self.prob_adjustment_table, \
                self.feature_extractor.get_feature_names())
        pickle.dump(state, file_handler)
        
    def load(self, file_handler):
        
        if self.model is not None or self.feature_selection is not None or self.default_prob_threshold is not None or \
                self.prob_adjustment_table is not None:

            log(('Warning: %s was already set with a model, feature selection, probability threshold or ' + \
                    'probability adjustment table; they will be overriden.') % self)
            
        self.model, self.feature_selection, self.default_prob_threshold, self.prob_adjustment_table, \
                feature_names = pickle.load(file_handler)
        assert self.feature_extractor.get_feature_names() == feature_names
        
    def _predict_proba_on_chunks(self, sample_chunks, thread_pool):

        n_samples = sum(map(len, sample_chunks))
        y_pred = np.array(n_samples * [np.nan])
        
        def _handle_chunk_features(chunk_data):
        
            chunk_start_index, chunk_X = chunk_data
                        
            if self.feature_selection is not None:
                chunk_X = self.feature_selection.transform(chunk_X)
                
            chunk_y_pred = predict_prob(self.model, chunk_X)
            y_pred[chunk_start_index:(chunk_start_index + len(chunk_y_pred))] = chunk_y_pred
        
        if thread_pool is None:
            self._run_feature_extraction_sequentially(sample_chunks, _handle_chunk_features)
        else:
            self._distribute_feature_extraction(sample_chunks, _handle_chunk_features, thread_pool)

        return y_pred
        
    def _run_feature_extraction_sequentially(self, sample_chunks, chunk_features_handler):
        
        chunk_start_index = 0
        
        for sample_chunk in sample_chunks:
            chunk_X = self.feature_extractor.create_features_as_matrix(sample_chunk)
            chunk_data = (chunk_start_index, chunk_X)
            chunk_features_handler(chunk_data)
            chunk_start_index += len(sample_chunk)
        
    def _distribute_feature_extraction(self, sample_chunks, chunk_features_handler, thread_pool):
        chunks_input_data = self._get_chunks_input_data(sample_chunks)
        distirbute_with_callback(thread_pool, invoke_feature_extractor_to_matrix, chunks_input_data, chunk_features_handler)
                    
    def _get_chunks_input_data(self, sample_chunks):
    
        chunk_start_index = 0
        chunks_input_data = []
        
        for sample_chunk in sample_chunks:
            chunks_input_data.append((self.feature_extractor_creator, self.feature_extractor_creator_args, self.feature_extractor_creator_kwargs, \
                    chunk_start_index, sample_chunk))
            chunk_start_index += len(sample_chunk)

        return chunks_input_data
        
    def _learn_prob_adjustment_table(self, X, y, bins):
        
        tune_mask = create_random_mask(len(X), len(X) // 2)
        X_tune, y_tune = X[tune_mask], y[tune_mask]
        X_train, y_train = X[~tune_mask], y[~tune_mask]

        self.model.fit(X_train, y_train)
        y_tune_pred = np.array(predict_prob(self.model, X_tune))
        
        self.prob_adjustment_table = []
        interval_length = 1 / bins
        
        for i in range(bins):
            min_value = i * interval_length
            max_value = (i + 1) * interval_length
            mask = (min_value <= y_tune_pred) & (y_tune_pred <= max_value)
            probability = y_tune[mask].mean()
            self.prob_adjustment_table.append(((min_value + max_value) / 2, probability))
            
class AsyncClassifier(object):
    
    def __init__(self, classifier_prediction_function, thread_pool, n_threads, chunk_size = 512):
        self.classifier_prediction_function = classifier_prediction_function
        self.thread_pool = thread_pool
        self.n_threads = n_threads
        self.chunk_size = chunk_size
        self.batch_size = n_threads * chunk_size
        self._active_submissions = []
        
    def submit_samples(self, samples, callback):
        
        self._active_submissions.append(_ClassificationSubmission(callback, list(samples)))
        n_pending_samples = self.n_pending_samples()
        log('Submitted %d new samples to classify. There are now %s pending samples.' % (len(samples), n_pending_samples))
        
        if n_pending_samples >= self.batch_size:
            self._process_sample_batch()
        
    def process_remaining_samples(self):
        while len(self._active_submissions) > 0:
            self._process_sample_batch()
        
    def n_pending_samples(self):
        return sum([len(submission.pending_samples) for submission in self._active_submissions])
        
    def close(self):
        self.thread_pool.close()
        
    def __enter__(self):
        pass
        
    def __exit__(self, ex_type, ex_value, tb):
        
        if ex_type is None:
            self.process_remaining_samples()
        
        self.close()
        
    def _process_sample_batch(self):
        batch_samples = self._get_next_batch()
        log('Classifying a batch of %d samples...' % len(batch_samples))
        batch_results = self.classifier_prediction_function(batch_samples, chunk_size = self.chunk_size, thread_pool = self.thread_pool)
        self._process_results(list(batch_results))
                
    def _get_next_batch(self):
    
        batch_samples = []
    
        for submission in self._active_submissions:
            if len(batch_samples) < self.batch_size:
                
                taken_samples = submission.pending_samples[:(self.batch_size - len(batch_samples))]
                batch_samples.extend(taken_samples)
                
                submission.pending_samples = submission.pending_samples[len(taken_samples):]
                submission.n_processing_samples += len(taken_samples)
            else:
                break
                
        return batch_samples
        
    def _process_results(self, results):
        
        n_submissions_to_remove = 0
        
        for submission in self._active_submissions:
            if len(results) > 0:
                
                submission_results = results[:submission.n_processing_samples]
                results = results[submission.n_processing_samples:]
                
                submission.n_processing_samples = 0
                submission.pending_results.extend(submission_results)
                
                if len(submission.pending_samples) == 0:
                    submission.callback(submission.pending_results)
                    n_submissions_to_remove += 1
            else:
                break
                
        self._active_submissions = self._active_submissions[n_submissions_to_remove:]
        
def predict_prob(model, X):
    return zip(*model.predict_proba(X))[1]
    
def find_prob_threshold(y_true, y_prob):
    
    '''
    Uses a binary search in order to find the threshold for y_prob such that the derived y_pred will result recall
    and specificity scores that are as close as possible.
    '''
    
    optional_thresholds = sorted(y_prob)
    threshold = None
    
    lower = 0
    upper = len(optional_thresholds) - 1
    
    while lower <= upper:
        
        i = (lower + upper) // 2
        threshold = optional_thresholds[i]
        
        y_pred = y_prob >= threshold
        recall = recall_score(y_true, y_pred)
        specificity = specificity_score(y_true, y_pred)
        
        if specificity < recall:
            # There are too few 0 predictions --> threshold is too low --> need to search for higher thresholds
            lower = i + 1
        else:
            upper = i - 1
            
    return threshold
        
class _ClassificationSubmission(object):
    def __init__(self, callback, samples):
        self.callback = callback
        self.pending_samples = samples
        self.n_processing_samples = 0
        self.pending_results = []