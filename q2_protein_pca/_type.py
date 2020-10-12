from q2_types.feature_data import FeatureData
from qiime2.core.type import SemanticType

RankedProteinAlignment = SemanticType('RankedProteinAlignment',
                                      variant_of=FeatureData.field['type'])

PositionMapping = SemanticType('PositionMapping',
                               variant_of=FeatureData.field['type'])
