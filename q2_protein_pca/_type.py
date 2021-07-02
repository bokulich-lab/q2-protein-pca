# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import FeatureData
from qiime2.core.type import SemanticType

RankedProteinAlignment = SemanticType('RankedProteinAlignment',
                                      variant_of=FeatureData.field['type'])

PositionMapping = SemanticType('PositionMapping',
                               variant_of=FeatureData.field['type'])
