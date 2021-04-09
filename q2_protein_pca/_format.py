# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv

from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model


class RankedProteinAlignmentFormat(model.TextFileFormat):
    HEADER = ["Sequence ID", "pos"]

    def _check_n_records(self, n=None):
        with self.open() as fh:
            reader = csv.reader(fh, delimiter=',')
            data_line_count = 0
            header = None
            header_len = None

            file_ = enumerate(reader) if n is None else zip(range(n), reader)

            for i, line in file_:
                # Tracks line number for error reporting
                i = i + 1

                if header is None:
                    header_first = line[0]
                    header_rest = set([x[:3] for x in line[1:]])
                    if header_first != self.HEADER[0]:
                        raise ValidationError(
                            '%s must be the first header value. The '
                            'first header values provided is: %s (on '
                            'line %s).' % (self.HEADER[0], header_first, i))
                    if header_rest != {self.HEADER[1]}:
                        raise ValidationError(
                            '"%s" must be the value of headers starting from '
                            'the second one. The first header values provided '
                            'is: "%s" (on line %s).' % (
                                self.HEADER[1], list(header_rest)[0], i))
                    header = [header_first]
                    header_len = len(line)
                else:
                    if len(line) != header_len:
                        raise ValidationError(
                            'Number of values on line %s are not the same as '
                            'number of header values. Found %s values '
                            '(%s), expected %s.' % (i, len(line), line,
                                                    header_len))

                    ranks = line[1:]
                    if not all([x.isdigit() for x in ranks]):
                        raise ValidationError(
                            'Some values on line %s are not numbers.' % i
                        )
                    if not all([(0 <= int(x) <= 23) for x in ranks]):
                        raise ValidationError(
                            'Some values on line %s are out of range.' % i
                        )

                    data_line_count += 1

            if data_line_count == 0:
                raise ValidationError('No taxonomy records found, only blank '
                                      'lines and/or a header row.')

    def _validate_(self, level):
        self._check_n_records(n={'min': 5, 'max': None}[level])


RankedProteinAlignmentDirectoryFormat = model.SingleFileDirectoryFormat(
    'RankedProteinAlignmentDirectoryFormat', 'ranked-protein-alignment.tsv',
    RankedProteinAlignmentFormat)


def _validate_record_min_len(cells, current_line_number, exp_len):
    if len(cells) < exp_len:
        raise ValidationError(
            "Expected data record to be TSV with ≥ {0} "
            "fields. Detected {1} fields at line {2}:\n\n{3!r}"
            .format(exp_len, len(cells), current_line_number, cells))


def _validate_file_not_empty(has_data):
    if not has_data:
        raise ValidationError(
            "There must be at least one data record present in the "
            "file in addition to the header line.")


class PositionMappingFormat(model.TextFileFormat):
    def _validate_(self, level):
        # n_records = {'min': 10, 'max': None}[level]
        with self.open():
            # # validate header
            # # for now we will not validate any information in the header.
            # line = fh.readline()
            #
            # # validate body
            # has_data = False
            # for line_number, line in enumerate(fh, start=2):
            #     cells = line.strip().split(',')
            #     _validate_record_min_len(cells, line_number, 2)
            #     for cell in cells[1:]:
            #         try:
            #             float(cell)
            #         except ValueError:
            #             raise ValidationError(
            #                 "Expected data to be comprised of float values. "
            #                 "Found non-float value {0} at line {1}"
            #                 .format(cell, line_number))
            #     has_data = True
            #     if n_records is not None and (line_number - 1) >= n_records:
            #         break
            has_data = True
            _validate_file_not_empty(has_data)


PositionMappingDirectoryFormat = model.SingleFileDirectoryFormat(
    'PositionMappingDirectoryFormat', 'position-mapping.csv',
    PositionMappingFormat)
