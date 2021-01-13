from .config import TRAINED_CLASSIFER_DUMP_FILE_PATH
from .variant_feature_extraction import FeatureExtractionSetup, get_snp_effect_feature_extractor
from .ml.classification import Classifier

def load_classifier(geneffect_setup, dump_file_path = TRAINED_CLASSIFER_DUMP_FILE_PATH):
    
    feature_extraction_setup = FeatureExtractionSetup(geneffect_setup)
    classifier = Classifier(get_snp_effect_feature_extractor, feature_extractor_creator_args = [feature_extraction_setup])
    
    with open(dump_file_path, 'rb') as f:
        classifier.load(f)
    
    return classifier