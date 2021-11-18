#include "grove.h"

    BranchCatalogue::BranchCatalogue(int levels)
    {
        if (levels > (1 << LEVEL_BITS))
        {
            logerr("Branch catalogue created with too many branch levels: %d. Max value is %d", levels, (1 << LEVEL_BITS) - 1);
            levels = (1 << LEVEL_BITS);
        }
        for (int i = 0; i < levels; i++)
        {
            branches.push_back(std::vector<PackedBranch>());
            occupied.push_back(BitVector());
            first_free.push_back(0);
        }
    }
    PackedBranch &BranchCatalogue::get(unsigned pos)
    {
        if (occupied[pos & ((1 << LEVEL_BITS) - 1)].get(pos >> LEVEL_BITS))
            return branches[pos & ((1 << LEVEL_BITS) - 1)][pos >> LEVEL_BITS];
        else
            return emptyBranch;
    }
    int BranchCatalogue::add(PackedBranch &b, int level)
    {
        if (level < 0 || level >= branches.size())
        {
            logerr("error adding branch level %d br size %d",level, branches.size());
            return -1;
        }
        int pos = first_free[level];
        if (pos == branches[level].size())
        {
            //increase vector size
            first_free[level]++;
            occupied[level].push_back(true);
            branches[level].push_back(b);
        }
        else
        {
            branches[level][pos] = b;
            occupied[level].set_unsafe(pos,true);
            first_free[level] = branches[level].size(); 
            for (int i=pos;i<branches[level].size();i++)
            {
                if (!occupied[level].get_unsafe(i))
                {
                    first_free[level] = i;
                    break;
                }
            }
        }
        return ((pos) << LEVEL_BITS) + level;
    }
    int BranchCatalogue::get_level_size(int level)
    {
        if (level < 0)
            level = 0;
        if (level >= branches.size())
            level = branches.size() - 1;
        int sz = 0;
        for (int i=0;i<occupied[level].size();i++)
        {
            sz += occupied[level].get_unsafe(i);
        }
        return sz;
    }
    void BranchCatalogue::remove(unsigned id)
    {
        int level = id & ((1 << LEVEL_BITS) - 1);
        int pos = id >> LEVEL_BITS;
        occupied[level].set(pos, false);
        first_free[level] = MIN(first_free[level],pos);
    }