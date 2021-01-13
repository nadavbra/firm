from collections import defaultdict, Counter

def setup_uniprot_tracks(geneffect_setup):

    all_uniprot_boolean_tracks = set()
    all_uniprot_categorical_tracks = defaultdict(set)
    
    for uniprot_record in geneffect_setup.uniprot_records.values():
        
        uniprot_record.boolean_tracks, uniprot_record.categorical_tracks = parse_tracks_from_uniprot_record(uniprot_record.raw_biopython_record)
        all_uniprot_boolean_tracks |= set(uniprot_record.boolean_tracks.keys())
        
        for track_name, track in uniprot_record.categorical_tracks.items():
            all_uniprot_categorical_tracks[track_name] = all_uniprot_categorical_tracks[track_name] | track.get_values()
        
    geneffect_setup.all_uniprot_boolean_tracks = list(sorted(all_uniprot_boolean_tracks))
    geneffect_setup.all_uniprot_categorical_tracks = list(sorted([(track_name, sorted(values)) for track_name, values in \
            all_uniprot_categorical_tracks.items()]))

def parse_tracks_from_uniprot_record(uniprot_biopython_record):
    
    boolean_tracks = defaultdict(lambda: Track(False))
    categorical_tracks = defaultdict(lambda: Track(None))
    
    for feature in uniprot_biopython_record.features:
        
        try:
            # Biopython converts the locations to 0-based indices. We need to restore it as 1-based.
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
        except TypeError:
            continue
        
        _add_to_boolean_tracks(boolean_tracks, _FEATURE_TO_BOOLEAN_TRACK_NAMES.get(feature.type, []), start, end)
        _add_to_categorical_tracks(categorical_tracks, _FEATURE_TO_CATEGORICAL_TRACKS.get(feature.type, []), start, end)
    
        if feature.type == 'modified residue':
            ptm_name = feature.qualifiers['description']
            _add_to_boolean_tracks(boolean_tracks, _get_ptm_track_names(ptm_name), start, end)
    
    return map(Track.simplify_track_dict, [boolean_tracks, categorical_tracks])
        
class Track(object):
        
    @staticmethod
    def simplify_track_dict(track_dict):
        
        for track in track_dict.values():
            track.simplify()
            
        return dict(track_dict)
    
    def __init__(self, default_value):
        self.default_value = default_value
        self._entries = []
        self._entries_by_value = defaultdict(list)
        
    def add_entry(self, start, end, value):
        entry = TrackEntry(start, end, value)
        self._entries.append(entry)
        self._entries_by_value[value].append(entry)
        
    def simplify(self):
        self._entries.sort(key = TrackEntry.start_getter)
        self._replace_with_simplified_entries()
        self._reset_entries_by_value()
        
    def calc_distance(self, location, total_length, value = True):
    
        value_entries = self._entries_by_value[value]
        
        if value == self.default_value:
            value_entries.extend(self._get_gap_entries(total_length))
    
        if len(value_entries) == 0:
            return None
        else:
            return min([entry.calc_distance(location) for entry in value_entries])
        
    def get_values(self):
        return set([entry.value for entry in self._entries])
        
    def count_values(self, total_length = None, start_bound = None, end_bound = None):
    
        def get_effective_start(start):
            if start_bound is None:
                return start
            else:
                return max(start, start_bound)
                
        def get_effective_end(end):
            if end_bound is None:
                return end
            else:
                return min(end, end_bound)
    
        counter = Counter()
        
        for entry in self._entries:
            
            effective_entry_length = get_effective_end(entry.end) - get_effective_start(entry.start) + 1
            
            if effective_entry_length > 0:
                counter[entry.value] += effective_entry_length
            
        if total_length is not None:
            
            effective_total_length = get_effective_end(total_length) - get_effective_start(1) + 1
            remained_length = effective_total_length - sum(counter.values())
            
            if remained_length > 0:
                counter[self.default_value] += remained_length
            
        return dict(counter)
        
    def _replace_with_simplified_entries(self):
    
        new_entries = []
    
        for entry in self._entries:
            
            new_entries.append(entry)
            
            if len(new_entries) >= 2:
                one_before_last_entry, last_entry = new_entries[-2:]
                new_entries = new_entries[:-2] + one_before_last_entry.merge_with(last_entry)
                
        self._entries = new_entries
        
    def _reset_entries_by_value(self):
        
        self._entries_by_value = defaultdict(list)
        
        for entry in self._entries:
            self._entries_by_value[entry.value].append(entry)
            
    def _get_gap_entries(self, total_length):
        
        last_entry_end = 0
            
        for entry in self._entries:
            if entry.start > last_entry_end + 1:
                yield TrackEntry(last_entry_end + 1, entry.start - 1, self.default_value)
                last_entry_end = entry.end
                
        if last_entry_end < total_length:
            yield TrackEntry(last_entry_end + 1, total_length, self.default_value)
    
    def __repr__(self):
        if len(self._entries) > 0:
            return ''.join([str(entry) for entry in self._entries])
        else:
            return '<empty>'
            
class TrackEntry(object):

    start_getter = lambda entry: entry.start

    def __init__(self, start, end, value):
        self.start = start
        self.end = end
        self.value = value
        self.length = self.end - self.start + 1
        self._validate_length()
        
    def calc_distance(self, location):
        if location < self.start:
            return self.start - location
        elif location <= self.end:
            return 0
        else:
            return location - self.end
        
    def touches(self, other):
        return max(self.start, other.start) <= min(self.end, other.end) + 1
        
    def merge_with(self, other):
        
        first_entry, second_entry = sorted([self, other], key = TrackEntry.start_getter)
    
        if first_entry.touches(second_entry):
            if first_entry.value == second_entry.value:
                return [TrackEntry(min(first_entry.start, second_entry.start), max(first_entry.end, second_entry.end), \
                        first_entry.value)]
            else:
                
                new_second_entry = TrackEntry(first_entry.end + 1, second_entry.end, second_entry.value)
                
                if new_second_entry.length > 0:
                    return [first_entry, new_second_entry]
                else:
                    return [first_entry]
        else:
            return [first_entry, second_entry]
    
    def __repr__(self):
        return '[%d-%d]{%s}' % (self.start, self.end, str(self.value))
        
    def _validate_length(self):
        if self.length <= 0:
            raise Exception('Invalid entry length: %d' % self.length)
            
def _add_to_boolean_tracks(tracks_dict, track_names, start, end, value = True):
    for track_name in track_names:
        tracks_dict[track_name].add_entry(start, end, value)
        
def _add_to_categorical_tracks(tracks_dict, tracks, start, end):
    for track_name, value in tracks:
        tracks_dict[track_name].add_entry(start, end, value)
        
def _get_ptm_track_names(raw_ptm_name):
    ptm_name = raw_ptm_name.split(';')[0]
    return _PTM_TO_BOOLEAN_TRACK_NAMES.get(ptm_name, []) + [_GLOBAL_PTM_TRACK]
    
_GLOBAL_PTM_TRACK = 'ptm'

_FEATURE_TO_BOOLEAN_TRACK_NAMES = {
    'active site': ['active-site'],
    'binding site': ['binding-site', 'generic-binding'],
    'DNA-binding region': ['binding-site', 'DNA-binding'],
    'lipid moiety-binding region': ['binding-site', 'lipid-binding'],
    'metal ion-binding site': ['binding-site', 'metal-binding'],
    'domain': ['domain', 'standard-domain'],
    'topological domain': ['domain', 'topological-domain'],
    'disulfide bond': ['disulfide', _GLOBAL_PTM_TRACK],
    'glycosylation site': ['glycosylation', _GLOBAL_PTM_TRACK],
}

_FEATURE_TO_CATEGORICAL_TRACKS = {
    'helix': [('ss', 'H')],
    'strand': [('ss', 'B')],
    'turn': [('ss', 'T')],
}

_PTM_TO_BOOLEAN_TRACK_NAMES = {

    'N-acetylalanine': ['acetyl', 'PTM-N'],
    'N-acetylaspartate': ['acetyl', 'PTM-N'],
    'N-acetylcysteine': ['acetyl', 'PTM-N'],
    'N-acetylglutamate': ['acetyl', 'PTM-N'],
    'N-acetylglycine': ['acetyl', 'PTM-N'],
    'N-acetylmethionine': ['acetyl', 'PTM-N'],
    'N-acetylserine': ['acetyl', 'PTM-N'],
    'N-acetylthreonine': ['acetyl', 'PTM-N'],
    'N-acetylvaline': ['acetyl', 'PTM-N'],
    'N6-acetyllysine': ['acetyl', 'PTM-N'],
    'O-acetylserine': ['acetyl', 'PTM-O'],
    
    'Diphthamide': ['amide'],
    'Glycine amide': ['amide'],
    'Isoleucine amide': ['amide'],
    'Leucine amide': ['amide'],
    'Methionine amide': ['amide'],
    'Phenylalanine amide': ['amide'],
    'Tyrosine amide': ['amide'],
    
    'O-AMP-threonine': ['AMP', 'PTM-O'],
    'O-AMP-tyrosine': ['AMP', 'PTM-O'],

    'N6-biotinyllysine': ['biotinyl', 'PTM-N'],
    
    'Deamidated asparagine': ['Deamidated'],
    'Deamidated glutamine': ['Deamidated'],

    'Beta-decarboxylated aspartate': [],
    '4-carboxyglutamate': ['carboxyl', 'PTM-4'],
    'N6-carboxylysine': ['carboxyl', 'PTM-N'],
    'Pyrrolidone carboxylic acid': ['carboxyl'],
    
    '(3R)-3-hydroxyasparagine': ['hydroxy', 'PTM-3', 'PTM-R'],
    '(3R)-3-hydroxyaspartate': ['hydroxy', 'PTM-3', 'PTM-R'],
    '(3S)-3-hydroxyasparagine': ['hydroxy', 'PTM-3', 'PTM-S'],
    '3-hydroxyasparagine': ['hydroxy', 'PTM-3'],
    '3-hydroxyproline': ['hydroxy', 'PTM-3'],
    '4-hydroxyproline': ['hydroxy', 'PTM-4'],
    '5-hydroxylysine': ['hydroxy', 'PTM-5'],
    'Hydroxyproline': ['hydroxy'],
    
    'S-8alpha-FAD cysteine': ['FAD', 'PTM-S'],
    'Tele-8alpha-FAD histidine': ['FAD'],
    
    'Nitrated tyrosine': ['nitro'],
    'S-nitrosocysteine': ['nitro', 'PTM-S'],
    
    'Phosphohistidine': ['phosphorylation'],
    'Phosphoserine': ['phosphorylation'],
    'Phosphothreonine': ['phosphorylation'],
    'Phosphotyrosine': ['phosphorylation'],
    
    'ADP-ribosylarginine': ['ADP-ribosyl'],
    'ADP-ribosylasparagine': ['ADP-ribosyl'],
    'ADP-ribosylcysteine': ['ADP-ribosyl'],
    
    'Asymmetric dimethylarginine': ['methyl', 'dimethyl'],
    'Cysteine methyl ester': ['methyl', 'monomethyl'],
    'Dimethylated arginine': ['methyl', 'dimethyl'],
    'N,N,N-trimethylalanine': ['methyl', 'trimethyl', 'PTM-N'],
    'N,N-dimethylproline': ['methyl', 'dimethyl', 'PTM-N'],
    'N4,N4-dimethylasparagine': ['methyl', 'dimethyl', 'PTM-N'],
    'N6,N6,N6-trimethyllysine': ['methyl', 'trimethyl', 'PTM-N'],
    'N6,N6-dimethyllysine': ['methyl', 'dimethyl', 'PTM-N'],
    'N6-methylated lysine': ['methyl', 'monomethyl', 'PTM-N'],
    'N6-methyllysine': ['methyl', 'monomethyl', 'PTM-N'],
    'Omega-N-methylarginine': ['methyl', 'monomethyl', 'PTM-N'],
    'Omega-N-methylated arginine': ['methyl', 'monomethyl', 'PTM-N'],
    'Pros-methylhistidine': ['methyl', 'monomethyl'],
    'S-methylcysteine': ['methyl', 'monomethyl', 'PTM-S'],
    'Symmetric dimethylarginine': ['methyl', 'dimethyl'],
    'Tele-methylhistidine': ['methyl', 'monomethyl'],
    
    'Cysteine persulfide': ['sulfide'],
    'Cysteine sulfenic acid (-SOH)': ['sulfide'],
    'Cysteine sulfinic acid (-SO2H)': ['sulfide'],
    'Methionine (R)-sulfoxide': ['sulfide'],
    'Methionine sulfoxide': ['sulfide'],
    'S-(dipyrrolylmethanemethyl)cysteine': ['PTM-S'],
    'Sulfoserine': ['sulfide'],
    'Sulfotyrosine': ['sulfide'],
    
    '1-thioglycine': ['PTM-1'],
    '2\',4\',5\'-topaquinone': ['PTM-2', 'PTM-4', 'PTM-5'],
    '3-oxoalanine (Cys)': ['PTM-3'],
    '5-glutamyl glycerylphosphorylethanolamine': ['glutamyl', 'PTM-5'],
    '5-glutamyl polyglutamate': ['glutamyl', 'PTM-5'],
    'Allysine': ['allysine'],
    'Blocked amino end (Ser)': [],
    'Citrulline': ['citrulline'],
    'Glycyl adenylate': [],
    'N6-(pyridoxal phosphate)lysine': ['pyridoxal', 'PTM-N'],
    'N6-(retinylidene)lysine': ['PTM-N'],
    'N6-glutaryllysine': ['glutary', 'PTM-N'],
    'N6-lipoyllysine': ['PTM-N'],
    'N6-malonyllysine': ['malonyl', 'PTM-N'],
    'N6-succinyllysine': ['succinyl', 'PTM-N'],
    'N-pyruvate 2-iminyl-valine': ['PTM-N'],
    'S-cysteinyl cysteine': ['PTM-S'],
    'Thyroxine': [],
    'Triiodothyronine': [],
}