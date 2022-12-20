# -*- coding: utf-8 -*-
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import pycountry
from .us_states import us_state_to_code, code_to_us_state


STATES = r'Alabama|Alaska|Arizona|Arkansas|California|Colorado|Connecticut|Delaware|Florida|Georgia|' \
         'Hawaii|Idaho|Illinois|Indiana|Iowa|Kansas|Kentucky|Louisiana|Maine|Maryland|Massachusetts|' \
         'Michigan|Minnesota|Mississippi|Missouri|Montana|Nebraska|Nevada|New[\s_]Hampshire|New[\s_]Jersey|' \
         'New[\s_]Mexico|New[\s_]York|North[\s_]Carolina|North[\s_]Dakota|Ohio|Oklahoma|Oregon|Pennsylvania|Rhode[\s_]Island|' \
         'South[\s_]Carolina|South[\s_]Dakota|Tennessee|Texas|Utah|Vermont|Virginia|Washington|West[\s_]Virginia|' \
         'Wisconsin|Wyoming'


class IAVSequence(object):
    """
    Represents an IAV genetic sequence (nucleotide or aa).
    The class parses a sequence's name; in particular, the tokens separated by '|'.
    All header tokens are then stored in self.tokens.

    Note that each sequence name must include a date and sequence identifier (e.g., 'A/swine/IA/12345/2015').
    The date must be in one of the following formats: %Y-%m-%d, %Y-%m, %m/%d/%Y, %m/%Y, %Y.
    For example, 2014-03-14, 2015-07, 03/05/2019, 01/2016, 2017 are all acceptable dates.
    """

    def __init__(self, full_name: str, seq: str, protein='HA', host='swine'):
        self.full_name = full_name
        self.seq = seq.lower()
        assert self.seq.count('-') != len(self.seq)  # Make sure there are informative sites.
        self.protein = protein
        self.host = host
        self.subtype = 'unknown'
        self.name = None
        self.date = None
        self.state = None
        self.country_code = None
        self.tokens = set()

        self._parse_sequence_header(full_name)
        usda_match = re.search(r'A0\d+', full_name)
        self.usda_id = None
        if usda_match:
            self.usda_id = usda_match.group()
        # try:
        #     if host == 'swine':
        #         self.strain_name, self.type, self.subtype, self.state, self.host,\
        #             pdm09_str, self.h1clade, date_str = full_name.split('|')
        #     else:
        #         self.strain_name, self.type, self.subtype, self.state, self.host, \
        #             pdm09_str, date_str = full_name.split('|')
        #     self.pdm09 = (pdm09_str == 'Y')
        #     self.host = self.host.lower()
        #     if self.host == 'swine':
        #         self.id = self.strain_name.split('/')[3]
        #     else:
        #         self.id = self.strain_name.split('/')[2]
        #     if len(date_str.split('/')) == 3:
        #         self.date = datetime.strptime(date_str, '%m/%d/%Y')
        #     elif len(date_str.split('/')) == 2:
        #         self.date = datetime.strptime(date_str, '%m/%Y')
        #     else:
        #         self.date = datetime.strptime(date_str, '%Y')
        #     self.parsed_name = True
        # except Exception as err:
        #     print('parsing error')
        #     print(err)
        #     self.parsed_name = False

    def _parse_sequence_header(self, full_name):
        for token in full_name.split('|'):
            if not token:
                continue
            self.tokens.add(token)
            if token.count('/') > 0 and self.matches(token, [r'.*[a-zA-Z].*']) and \
                    self.matches(token, [r'[\w\d/_\-\(\)]+']):
                # Name: contains '/', contains a letter, consists of letter/digit/'/'/_/- combination.
                self.name = token
            elif self.matches(token, [r'[A-Z]{3}']) and pycountry.countries.get(alpha_3=token):
                # Country code.
                self.country_code = token
            elif self.matches(token, [r'[\d\-/]{4,}']) and not self.matches(token, [r'\d{5,}']):
                # Date: consists of digits and '/' or '-' symbols.
                try:
                    self.date_str = token
                    if token.count('/') == 2:
                        self.date = datetime.strptime(token, '%m/%d/%Y')
                    elif token.count('/') == 1:
                        self.date = datetime.strptime(token, '%m/%Y')
                    elif token.count('-') == 2:
                        self.date = datetime.strptime(token, '%Y-%m-%d')
                    elif token.count('-') == 1:
                        self.date = datetime.strptime(token, '%Y-%m')
                    else:
                        self.date = datetime.strptime(token, '%Y')
                except Exception as e:
                    print(e)
                    raise Exception('Cannot parse date %s in header %s' % (token, full_name))
            elif self.matches(token.lower(), ['swine|human']):
                # Host.
                self.host = token.lower()
            elif self.matches(token.upper(), [r'H\dN\d|H\d|N\d']):
                # Subtype.
                self.subtype = token.upper()
            elif len(token) == 1 and self.matches(token.upper(), [r'A|B|C|D']):
                # Type.
                self.type = token.upper()
            elif token in us_state_to_code.keys():
                self.state = token
            elif token in code_to_us_state.keys():
                self.state = code_to_us_state[token]
        if not self.name or not self.date:
            raise Exception('The strain name should contain a unique identifier and a date: ' + full_name)

    def matches(self, token, regexes: [str]):
        for regex in regexes:
            if re.fullmatch(regex, token, flags=re.ASCII):
                return True
        return False

    def is_swine(self) -> bool:
        if self.host:
            return self.host == 'swine'
        else:
            return self.full_name.lower().count('swine') > 0

    def is_human(self) -> bool:
        if self.host:
            return self.host.lower() == 'human'
        else:
            return not self.is_swine()

    def seq_start(self):
        return re.search(r'[^-]', str(self.seq)).start()

    def seq_end(self):
        return len(self.seq) - re.search(r'[^-]', str(self.seq[::-1])).start() - 1

    def to_bio(self) -> SeqRecord:
        """
        Represents IAVSequence as an instance of Bio.SeqRecord (from biopython).
        """
        if isinstance(self.seq, Seq):
            return SeqRecord(self.seq, id=self.full_name, description='')
        else:
            return SeqRecord(Seq(str(self.seq)), id=self.full_name, description='')

    def contains_token(self, token: str) -> bool:
        return token in self.tokens

    def in_header(self, substr: str) -> bool:
        return self.full_name.count(substr) > 0

    def get_swine_name_id(self):
        name_tokens = self.name.split('/')
        if len(name_tokens) == 5:
            return name_tokens[3]
        else:
            return None

    def __repr__(self):
        # strain_str = '%s\n Name: %s\n Date: %s\n Subtype: %s' %\
        #              (self.full_name, self.name, self.date, self.subtype)
        return self.full_name
