#!/usr/bin/env python3

"""
Add functionality for EC record parsing to AMON
"""


from KEGG_parser.parsers import *
from KEGG_parser.downloader import get_from_kegg_flat_file, get_kegg_record_dict
from functools import partial
import re
import sys


### functions copied and modded from KEGG_parser
def split_entry(current_dict, current_entry_name, current_entry_data, field=0):
    """replaces KEGG_parser split_entry by adding field arg because EC numbers 
       are all prefixed with "EC ", so we want to return second field
    """
    current_dict[current_entry_name] = current_entry_data.split()[field]
    return(current_dict)

def split_reaction(current_dict, current_entry_name, current_entry_data):
    """return ALL_REAC field in ec record as list
    Skip over reactions starting with "(other)" and remove symbols > and ;
    because comparable KO record lists these reactions 
    without any symbols, and no "(other)" ones 
    """
    if current_entry_name not in current_dict:
         current_dict[current_entry_name] = list()
    if current_entry_data.startswith("(other)"):
        return(current_dict)
    rxns = re.sub('[^\w ]+', '', current_entry_data).split()
    current_dict[current_entry_name].extend(rxns)
    return(current_dict)

#: how to parse each EC record
PARSE_EC_BY_FIELD = {
    'ENTRY': partial(split_entry, field=1), 'NAME': split_name_by_comma, 'ALL_REAC': split_reaction,
    'PATHWAY': split_and_append, 'MODULE': split_and_append, 'DISEASE': split_and_append,
    'CLASS': add_class, ##'DBLINKS': add_nested_dict, 'GENES': add_nested_dict
}


def parse_ec(ec_raw_record):
    """parse ec record; mod from `KEGG_parser.parsers.parse_ko`"""
    ec_dict = dict()
    past_entry = None
    for line in ec_raw_record.strip().split('\n'):
        current_entry_name = line[:12].strip()

        if current_entry_name == '':
            current_entry_name = past_entry

        current_entry_data = line[12:].strip()

        if current_entry_name != '':
            if current_entry_name in PARSE_EC_BY_FIELD:
                ec_dict = PARSE_EC_BY_FIELD[current_entry_name](ec_dict, current_entry_name, current_entry_data)
        past_entry = current_entry_name
    return ec_dict

def get_rns_from_ec_dict(ec_dict: dict):
    """get reactions from dictionary output by
       `KEGG_parser.downloader.get_kegg_record_dict`
    """
    list_of_ecs = ec_dict.keys()
    rns = set()
    for ec in list_of_ecs:
        try:
            ec_record = ec_dict[ec]
            if 'ALL_REAC' in ec_record:
                rns.update(ec_record['ALL_REAC'])
        except KeyError:
            pass
    return rns

