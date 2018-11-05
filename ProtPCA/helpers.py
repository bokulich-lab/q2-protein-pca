import os
import sys
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez


def set_logger(verbosity, log_to_stdout=False, log_file=None,
               log_format='%(asctime)s %(levelname)s\t%(message)s'):
    """ Sets a logger for the micro service
    Args:
        verbosity (int): Level of verbosity for the logger (1,2,3)
        log_to_stdout (boolean, optional): Write to standard system output
        log_file (:obj:`str`, optional): Name of the file to write the logs
        log_format (:obj:`str`, optional): Format of the logger output
    Returns:
        (:obj:logging.Logger) Logger object
    """
    verbosity_levels = {1: logging.ERROR, 2: logging.INFO, 3: logging.DEBUG}
    log_level = verbosity_levels[verbosity]
    if log_to_stdout:
        h = logging.StreamHandler(stream=sys.stdout)
    else:
        if not log_file:
            # Use module name
            log_file = os.path.splitext(os.path.split(__file__)[-1])[0]+'.log'
        h = logging.FileHandler(log_file)
    lf = logging.Formatter(log_format)
    h.setFormatter(lf)
    logger = logging.getLogger()
    logger.addHandler(h)
    logger.setLevel(log_level)
    return logger


class fasta:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence


def read_fasta(file):
    items = []
    index = 0
    for line in file:
        if line.startswith(">"):
            if index >= 1:
                items.append(aninstance)
            index += 1
            name = line[:-1]
            seq = ''
            aninstance = fasta(name, seq)
        else:
            seq += line[:-1]
            aninstance = fasta(name, seq)
    items.append(aninstance)
    return items


def split_string(text):
    return [text[i] for i in range(len(text))]


def find_unique_aas(df):
    list_of_unique_aas = []
    columns = df.columns.values.tolist()
    for i in range(len(df.columns)):
        list_of_unique_aas.append(df[columns[i]].unique())
    return list_of_unique_aas


def calculate_covariance_matrix(df):
    n = df.shape[0]  # number of rows in the dataframe (no. of sequences)
    m3 = np.ones((n, n))  # np array filled with ones that corresponds to the term 11'

    dev_scores = df.values - (m3.dot(df.values)) / n
    a_a = np.dot(np.transpose(dev_scores), dev_scores)
    v_co = a_a / n

    return v_co  # returns a variance-covariance matrix


def save_pos_to_file(alist, location, title):
    f = open(location, "w")
    f.write(title + "\n")
    for i in range(len(alist)):
        f.write(str(alist[i]) + "\n")
    f.close()


def euclidean_distance(x, y):
    return np.sqrt(x**2 + y**2)


def make_label(number, trail_digits=4):
    return '{:>0{trailing}}'.format(number, trailing=trail_digits)


# def get_positions(list_with_positions, list_with_values):
#     l = []
#     for i in range(len(list_with_positions)):
#         l.append([list_with_positions[i], list_with_values[list_with_positions[i]]])
#     return l


def locate_positions(x_coords, y_coords, conserved_coff=0.1, variable_coff=0.4):
    e_distances = [euclidean_distance(x, y) for (x, y) in zip(x_coords, y_coords)]
    max_e = max(e_distances)
    conserved_positions, variable_positions, other_positions= [], [], []

    for i, _ in enumerate(x_coords):
        x, y = x_coords[i], y_coords[i]
        if e_distances[i] < conserved_coff * max_e:
            conserved_positions.append([i, x, y])
        elif e_distances[i] > variable_coff * max_e:
            variable_positions.append([i, x, y])
        else:
            other_positions.append([i, x, y])
    return conserved_positions, variable_positions, other_positions


def positions_to_color(list_of_variable_positions, aligned_df):
    # remove gaps but preserve original indexing
    raw_df = aligned_df[aligned_df != '-']
    # get aa positions from the original sequence (no gaps) numbered using indices from the alignment
    positions = [int(x[3:]) for x in raw_df.index]
    variable_gap_cnt = 0
    variable_pos_no_gaps = []
    for j, var_pos in enumerate(list_of_variable_positions):
        m, m_val = var_pos[0], var_pos[1]
        pos2check = f'pos{m}'
        if aligned_df[pos2check] == "-":
            variable_gap_cnt += 1
        else:
            n_pos = positions.index(m)
            n_val = raw_df.loc[pos2check]
            variable_pos_no_gaps.append([m + 1, n_pos + 1, n_val, m_val, np.absolute(m_val)])
    return pd.DataFrame(variable_pos_no_gaps, columns=['PosAln', 'PosRaw', 'AA', "PC score", "Abs PC score"])


def fetch_record_ncbi(gi_number, email, db='protein', rettype='gb', retmode='text'):

    Entrez.email = email  # Always tell NCBI who you are
    handle = Entrez.efetch(db=db, id=gi_number, rettype=rettype, retmode=retmode)
    return handle