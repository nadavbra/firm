from __future__ import absolute_import, division, print_function

import traceback
from StringIO import StringIO
from collections import Counter
from itertools import chain

import numpy as np
import pandas as pd

from .util import log, IDENTITY_FUNCTION

def invoke_feature_extractor_to_matrix(data, verbose = True):
    try:
    
        feature_extractor_creator, feature_extractor_creator_args, feature_extractor_creator_kwargs, meta_data, samples = data
        feature_extractor = feature_extractor_creator(*feature_extractor_creator_args, **feature_extractor_creator_kwargs)
    
        if verbose:
            log('Extracting features from [%s]...' % meta_data)
    
        return meta_data, feature_extractor.create_features_as_matrix(samples)
    except:
        # Since this method is run through the multiprocessing.Pool interface, only the exception message is shown, and not
        # the stacktrace. 
        exception_str = StringIO()
        traceback.print_exc(file = exception_str)
        exception_str.seek(0)
        raise Exception(exception_str.read())

class FeatureExtractor(object):
        
    def n_features(self):
        '''
        @return: The number of features this extractor creates (int).
        '''
        raise NotImplementedError()
        
    def get_feature_names(self):
        '''
        @return: The names of the features this extractor creates (list of strings).
        '''
        raise NotImplementedError()
        
    def create_features_as_dict(self, sample):
        '''
        Creates features from the given sample object.
        @return: A dict of the following format:
        {
            'feature_name1': feature_value,
            'feature_name2': feature_value,
            ...
        }
        '''
        raise NotImplementedError()
        
    def create_features_into_array(self, feature_array, sample):
        '''
        Creates features from the given sample object, and puts them, by their order (as described by get_feature_names)
        into the given array (which must be of the exact size of n_features).
        '''
        raise NotImplementedError()
        
    def create_features_as_array(self, sample):
        '''
        Creates features from the given sample object.
        @return: The created features, as a float numpy array.
        '''
        feature_array = np.empty(self.n_features(), dtype = np.float64)
        self.create_features_into_array(feature_array, sample)
        return feature_array
        
    def create_features_as_matrix(self, samples, show_progress_bar = False):
        
        '''
        Creates featres for all the given sample objects.
        @return: The created features, as a float numpy matrix (shape: n_samples X n_features).
        '''
        
        if show_progress_bar:
            
            from IPython.display import display
            from ipywidgets import FloatProgress
            
            progress_bar = FloatProgress(min = 0, max = len(samples) - 1)
            display(progress_bar)

        feature_matrix = np.empty((len(samples), self.n_features()), dtype = np.float64)
        
        for i, sample in enumerate(samples):
            
            self.create_features_into_array(feature_matrix[i, :], sample)
            
            if show_progress_bar:
                progress_bar.value = i
        
        return feature_matrix
        
    def create_features_as_dataframe(self, samples, show_progress_bar = False):
        '''
        Creates featres for all the given sample objects.
        @return: The created features, as a Pandas dataframe with float values (row for each sample, column for each feature).
        '''
        return pd.DataFrame(self.create_features_as_matrix(samples, show_progress_bar = show_progress_bar), \
                columns = self.get_feature_names())

class CombinerFeatureExtractor(FeatureExtractor):

    def __init__(self, feature_extractors):
        
        self.feature_extractors = feature_extractors
        
        self._n_features_per_extractor = [feature_extractor.n_features() for feature_extractor in self.feature_extractors]
        self._feature_names_per_extractor = [feature_extractor.get_feature_names() for feature_extractor in self.feature_extractors]
        assert map(len, self._feature_names_per_extractor) == self._n_features_per_extractor
        
        self._n_features = sum(self._n_features_per_extractor)
        self._feature_names = [feature_name for feature_names in self._feature_names_per_extractor for \
                feature_name in feature_names]
        assert len(self._feature_names) == self._n_features
        
    def n_features(self):
        return self._n_features
        
    def get_feature_names(self):
        return self._feature_names
        
    def create_features_as_dict(self, sample):
        
        feature_dict = {}
        
        for feature_extractor in self.feature_extractors:
            feature_dict.update(feature_extractor.create_features_as_dict(sample))
        
        assert len(feature_dict) == self._n_features
        assert set(feature_dict.keys()) == set(self._feature_names)
        return feature_dict
        
    def create_features_into_array(self, feature_array, sample):
        
        index = 0
        
        for feature_extractor, n_features in zip(self.feature_extractors, self._n_features_per_extractor):
            feature_extractor.create_features_into_array(feature_array[index:(index + n_features)], sample)
            index += n_features

class SimpleFeatureExtractor(FeatureExtractor):

    def __init__(self, feature_names_or_extractors, extraction_function):
        
        '''
        @param feature_names_or_extractors (list): A list where each element is either a feature name (str) or feature extractor
        (FeatureExtractor). If an element is a string, the corresponding feature (outputed by extraction_function) will be
        considered a simple feature and taken as it is, with the provided feature name. If it's a feature extractor, the
        feature extractor will be used to extract the final features from the corresponding feature value outputed by 
        extraction_function.
        @param extraction_function (function): A function that takes a sample object and returns a sequence of features,
        corresponding to feature_names_or_extractors (if there's only one feature, can just return it).
        '''
        
        self._extraction_function = extraction_function
        self._feature_names_or_extractors = feature_names_or_extractors
        self._final_feature_names = list(chain(*[[feature] if isinstance(feature, str) else feature.get_feature_names() for \
                feature in feature_names_or_extractors]))
        self._set_feature_transformations()
        
    def n_features(self):
        return len(self._final_feature_names)
        
    def get_feature_names(self):
        return self._final_feature_names
        
    def create_features_as_dict(self, sample):
    
        feature_dict = {}
        
        for feature, value in zip(self._feature_names_or_extractors, self._extract_top_level_features(sample)):
            if isinstance(feature, str):
                feature_dict[feature] = value
            else:
                feature_dict.update(feature.create_features_as_dict(value))
    
        return feature_dict
        
    def create_features_into_array(self, feature_array, sample):
    
        top_level_feature_values = self._extract_top_level_features(sample)
        top_level_index = 0
        output_index = 0
        
        for feature_transformation in self._feature_transformations:
            if isinstance(feature_transformation, int):
                n_features = feature_transformation
                feature_array[output_index:(output_index + n_features)] = \
                        map(float, top_level_feature_values[top_level_index:(top_level_index + n_features)])
                top_level_index += n_features
                output_index += n_features
            else:
                inner_extractor, output_size = feature_transformation
                relevant_feature_array = feature_array[output_index:(output_index + output_size)]
                inner_extractor.create_features_into_array(relevant_feature_array, top_level_feature_values[top_level_index])
                top_level_index += 1
                output_index += output_size
        
    def _extract_top_level_features(self, sample):
        
        extracted = self._extraction_function(sample)
        
        try:
            return list(extracted)
        except TypeError:
            return [extracted]
        
    def _set_feature_transformations(self):
        
        consequent_simple_features = 0
        self._feature_transformations = []
        
        for feature in self._feature_names_or_extractors:
            if isinstance(feature, str):
                consequent_simple_features += 1
            else:
                
                if consequent_simple_features > 0:
                    self._feature_transformations.append(consequent_simple_features)
                    consequent_simple_features = 0
                
                self._feature_transformations.append((feature, feature.n_features()))
        
        if consequent_simple_features > 0:
            self._feature_transformations.append(consequent_simple_features)
            
class OneHotEncodingFeatureExtractor(FeatureExtractor):

    def __init__(self, feature_name_prefix, optional_values, value_extraction_function = IDENTITY_FUNCTION):
        
        self.feature_name_prefix = feature_name_prefix
        self.optional_values = optional_values
        self.value_extraction_function = value_extraction_function
        
        self.value_to_index = {value: i for i, value in enumerate(optional_values)}
        self.final_feature_names = [feature_name_prefix + str(value) for value in optional_values]
        
    def n_features(self):
        return len(self.optional_values)
        
    def get_feature_names(self):
        return self.final_feature_names
        
    def create_features_as_dict(self, sample):
        sample_value = self.value_extraction_function(sample)
        return {self.feature_name_prefix + str(value): (1 if value == sample_value else 0) for value in self.optional_values}
        
    def create_features_into_array(self, feature_array, sample):
        
        sample_value = self.value_extraction_function(sample)
        feature_array[:] = 0
        
        if sample_value in self.value_to_index:
            feature_array[self.value_to_index[sample_value]] = 1
        
class CounterFeatureExtractor(FeatureExtractor):

    def __init__(self, feature_name_prefix, optional_values, value_extraction_function = IDENTITY_FUNCTION):
        
        self.feature_name_prefix = feature_name_prefix
        self.optional_values = optional_values
        self.value_extraction_function = value_extraction_function
        
        self.value_to_index = {value: i for i, value in enumerate(optional_values)}
        self.final_feature_names = [feature_name_prefix + str(value) for value in optional_values]
        
    def n_features(self):
        return len(self.optional_values)
        
    def get_feature_names(self):
        return self.final_feature_names
        
    def create_features_as_dict(self, sample):
        sample_values = self.value_extraction_function(sample)
        return {self.feature_name_prefix + str(value): count for value, count in Counter(sample_values).items()}
        
    def create_features_into_array(self, feature_array, sample):
        
        sample_values = self.value_extraction_function(sample)
        feature_array[:] = 0
        
        for value in sample_values:
            if value in self.value_to_index:
                feature_array[self.value_to_index[value]] += 1