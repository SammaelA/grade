#include "function_stat.h"
#include "common_utils/distribution.h"
#include "common_utils/utility.h"

void my_opt::get_function_stat(std::vector<float> &param_list, const OptFunction &my_opt_f, FunctionStat &stat, bool reload)
{
    int batch_size = 512;
    int batches = 300;
    stat = FunctionStat(param_list.size());
    stat.name = my_opt_f.name;
    stat.version = my_opt_f.version;
    
    std::vector<std::vector<float>> params(batch_size, param_list);

    BlkManager man;
    Block b1;
    if (!reload)
        man.load_block_from_file("function_values_stat.blk", b1);
    Block *b2 = b1.get_block(stat.name);
    if (b2 && stat.version == b2->get_int("version",-1))
        stat.save_load_blk(*b2, false);
    else
        b2 = new Block();
    for (int b=0;b<batches;b++)
    {
        for (auto &par : params)
            for (auto &p : par)
                p = urand(0,1);
        auto res = my_opt_f.f(params);
        stat.tries += params.size();
        for (int j=0;j<res.size();j++)
        {
            float r = res[j];
            int q = 0;
            if (r < 0.001)
                q = 0;
            else
                q = CLAMP((int)(r*stat.Q_NUM) + 1,1,stat.Q_NUM);
            stat.quantiles[q]++;
            float mark = 0;
            for (int k=0;k<param_list.size();k++)
            {
                int bucket = CLAMP(params[j][k]*stat.Q_NUM,0,stat.Q_NUM-1);
                stat.by_variable_values[k][bucket][q]++;
                mark += stat.marks[k][bucket];
            }
            if (b > batches/2)
            {
                int mq = CLAMP(mark + 10,0,19);
                stat.marks_q[mq][0] += r > 0.001;
                stat.marks_q[mq][1] ++;
            }
        }

        float global_nzv = 1 - stat.quantiles[0]/(float)stat.tries;
        int i=0;
        for (auto &var : stat.by_variable_values)
        {
            for (int j=0;j<stat.Q_NUM;j++)
            {
                int z_cnt = var[j][0];
                int cnt = 0;
                for (int k=0;k<stat.Q_NUM+1;k++)
                    cnt += var[j][k];
                float nzv = 1 - (float)(z_cnt)/cnt;
                float mark = 0;
                if (nzv < 1e-6)
                    mark = -9.9;
                else if (nzv < global_nzv)
                    mark = 1-global_nzv/nzv;
                else
                    mark = nzv/global_nzv-1;
                stat.marks[i][j] = mark;
            }
            i++;
        }
        //stat.print();
    }

    b2->set_int("version", stat.version);
    stat.save_load_blk(*b2, true);
    b1.set_block(stat.name, b2);
    man.save_block_to_file("function_values_stat.blk", b1);
}
my_opt::FunctionStat::FunctionStat(int variables_cnt)
{
    if (variables_cnt == 0)
        return;
    by_variable_values = std::vector<int[Q_NUM][Q_NUM+1]>(variables_cnt);
    marks = std::vector<float[Q_NUM]>(variables_cnt);
    for (int i=0;i<2*Q_NUM;i++)
    {
        marks_q[i][0] = 0;
        marks_q[i][1] = 0;
    }
    for (int i=0;i<Q_NUM+1;i++)
        quantiles[i] = 0;
}

void my_opt::FunctionStat::save_load_blk(Block &b, bool save)
{
    logerr("save load blk");
    //int version = b.get_int("version",0);
    if (save)
    {
        b.set_int("tries", tries);

        std::vector<int> qn;
        for (int i=0;i<Q_NUM+1;i++)
            qn.push_back(quantiles[i]);
        b.set_arr("quantiles", qn);

        Block *bl = new Block();
        for (auto &v : by_variable_values)
        {
            std::vector<int> bvv; 
            for (int i=0;i<Q_NUM;i++)
                for (int j=0;j<Q_NUM+1;j++)
                    bvv.push_back(v[i][j]);
            bl->add_arr("values", bvv);
        }
        b.set_block("by_variable_values", bl);
        delete bl;

        Block *m_bl = new Block();
        for (auto &v : marks)
        {
            std::vector<float> mvv;
            for (int i=0;i<Q_NUM;i++)
                mvv.push_back(v[i]);
            m_bl->add_arr("values", mvv);
        }
        b.set_block("marks", m_bl);
        delete m_bl;
        
        std::vector<int> marks_q_arr; 
        for (auto i = 0;i<2*Q_NUM;i++)
        {
            marks_q_arr.push_back(marks_q[i][0]);
            marks_q_arr.push_back(marks_q[i][1]);
        }
        b.set_arr("marks_q", marks_q_arr);
    }
    else
    {
        tries = b.get_int("tries");

        std::vector<int> qn;
        b.get_arr("quantiles", qn);        
        for (int i=0;i<MIN(qn.size(), Q_NUM+1);i++)
            quantiles[i] = qn[i];
        for (int i=qn.size(); i<Q_NUM+1;i++)
            quantiles[i] = 0;
        
        Block *bl = b.get_block("by_variable_values");
        if (bl)
        {
            for (int i=0;i<bl->size();i++)
            {
                std::vector<int> bvv; 
                bl->get_arr(i, bvv);
                if (bvv.size() != Q_NUM*(Q_NUM + 1))
                {
                    logerr("variable %d by_variable_values corrupted size %d", i, bvv.size());
                }
                else
                {
                    for (int j=0;j<bvv.size();j++)
                    {
                        by_variable_values[i][j / (Q_NUM + 1)][j % (Q_NUM + 1)] = bvv[j];
                    }
                }
            }
        }

        Block *m_bl = b.get_block("marks");
        if (m_bl)
        {
            for (int i=0;i<m_bl->size();i++)
            {
                std::vector<float> bvv; 
                m_bl->get_arr(i, bvv);
                if (bvv.size() != Q_NUM)
                {
                    logerr("variable %d marks corrupted", i);
                }
                else
                {
                    for (int j=0;j<bvv.size();j++)
                    {
                        marks[i][j] = bvv[j];
                    }
                }
            }
        }

        std::vector<int> marks_q_arr; 
        b.get_arr("marks_q", marks_q_arr);
        if (marks_q_arr.size() != 4*Q_NUM)
        {
            logerr("marks_q corrupted");
        }
        else
        {
            for (int i=0;i<marks_q_arr.size();i++)
            {
                marks_q[i/2][i%2] = marks_q_arr[i];
            }
        }
    }
}

void my_opt::FunctionStat::print()
{
    debug("Function \"%s\" stat:\n", name.c_str());
    debug("Value quantiles:\n");
    int aggregate = 1;
    for (int i=0;i<Q_NUM/aggregate + 1;i++)
    {
        int cnt = 0;
        if (i == 0)
            cnt = quantiles[0];
        else
        {
            for (int j=aggregate*(i-1)+1;j<=aggregate*i;j++)
               cnt += quantiles[j]; 
        }
        int d_cnt = ceil(log2(1 + cnt));
        if (i == 0)
            debug("%.2f      ", i/(float)Q_NUM);
        else
            debug("%.2f-%.2f ", (aggregate*(i-1))/(float)Q_NUM, aggregate*i/(float)Q_NUM);
        for (int j = 0;j<d_cnt;j++)
        {
            debug("#");
        }
        debug(" %d (%.1f%)\n",cnt,100*cnt/(float)tries);
    }

    debug("By variable values:\n");
    int i=0;
    for (auto &var : by_variable_values)
    {
        debug("x_%d: ",i);
        for (int j=0;j<Q_NUM;j++)
        {
            float mark = marks[i][j];
            if (mark > 0)
                debug("+%.1f ", mark);
            else
                debug("%.1f ", mark);
        }
        debugnl();
        i++;
    }

    debug("marks quantiles:\n");
    for (int i=0;i<20;i++)
    {
        if (marks_q[i][1] > 0)
            debug("%.1f -- %d %.2f\n", (float)(i-10), marks_q[i][1], (float)marks_q[i][0]/marks_q[i][1]);
    }
}