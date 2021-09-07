#include "clustering.h"

bool ClusterData::is_valid()
{
    return (base && !ACDA.originals.empty() &&
            ACDA.originals.size() == ACDA.rotations.size() &&
            ACDA.originals.size() == IDA.centers_par.size() &&
            ACDA.originals.size() == IDA.centers_self.size() &&
            ACDA.originals.size() == IDA.transforms.size() &&
            ACDA.originals.size() == IDA.type_ids.size() &&
            ACDA.originals.size() == IDA.tree_ids.size());
}