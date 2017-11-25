from __future__ import absolute_import, division, print_function

from scipy.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests

from sklearn.base import BaseEstimator
from sklearn.feature_selection.base import SelectorMixin
from sklearn.ensemble import RandomForestClassifier

class RandomForestClassifierWithCoef(RandomForestClassifier):
    '''
    From http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn
    '''
    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_
        
class CorrelationFeatureSelection(BaseEstimator, SelectorMixin):
        
    def __init__(self, threshold = 0.05):
        self.threshold = threshold

    def _get_support_mask(self):
        return self.qvals <= self.threshold

    def fit(self, X, y):
        
        self.rhos = []
        self.pvals = []
        
        n, d = X.shape
        
        for i in xrange(d):
            x = X[:, i]
            rho, pval = spearmanr(x, y)
            self.rhos += [rho]
            self.pvals += [pval]
            
        _, self.qvals, _, _ = multipletests(self.pvals, method = 'fdr-bh')
        
        self.rhos = np.array(self.rhos)
        self.pvals = np.array(self.pvals)
        self.qvals = np.array(self.qvals)
