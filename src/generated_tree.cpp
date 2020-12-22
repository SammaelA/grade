
#include "generated_tree.h"
#include "tinyEngine/utility.h"
#include "branch_clusterization.h"
#include <math.h>
#include <algorithm>
#define PI 3.14159265f
int seg_count = 0;
int iter = 0;
float sum_feed[10];
float count_feed[10];
bool dice(float x, float base)
{
    float f = (float)rand()/RAND_MAX;
    return f < (x/base);
}
bool dice_min(float x, float min)
{
    if (x<1e-8)
        x = 1e-8;
    return dice(min,x);
}
glm::vec3 TreeGenerator::rand_dir()
{
    glm::vec3 dir;
    dir.x = (2.0*rand())/RAND_MAX - 1.0;
    dir.y = (2.0*rand())/RAND_MAX - 1.0;
    dir.z = (2.0*rand())/RAND_MAX - 1.0;
    return glm::normalize(dir);
}
glm::vec3 rand_planar_dir()
{
    glm::vec3 dir;
    dir.x = (2.0*rand())/RAND_MAX - 1.0;
    dir.y = (2.0*rand())/RAND_MAX - 1.0;
    dir.z = (2.0*rand())/RAND_MAX - 1.0;
    return glm::normalize(dir);
}
void TreeGenerator::new_joint(Branch *b, Joint &j)
{
    if (b->level == curParams.max_depth() -1)
    {
        j.type = j.LEAF;
        Leaf *l = curTree.leaves->new_leaf();
        l->pos = j.pos;
        glm::vec3 rd1 = rand_dir();
        glm::vec3 rd2 = rand_dir();

        
        glm::vec3 a = j.pos + curParams.leaf_size_mult()*rd1 + 0.5f*curParams.leaf_size_mult()*rd2;
        glm::vec3 b = j.pos + 0.5f*curParams.leaf_size_mult()*rd2;
        glm::vec3 c = j.pos - 0.5f*curParams.leaf_size_mult()*rd2;
        glm::vec3 d = j.pos + curParams.leaf_size_mult()*rd1 - 0.5f*curParams.leaf_size_mult()*rd2;
        ////fprintf(stderr,"(%f %f %f)(%f %f %f)(%f %f %f)\n",a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z);
        l->edges.push_back(a);
        l->edges.push_back(b);
        l->edges.push_back(c);
        l->edges.push_back(d);

        j.leaf = l;
    }
    else if (b->joints.empty())
    {
        j.type = j.MIDDLE;
    }
    else 
    {
        float b_ch = powf((b->base_seg_n +(float)b->joints.size()/(b->base_seg_n + b->max_seg_count)), 
                          curParams.branching_power());
        b_ch = curParams.min_branching_chance() + b_ch*(curParams.max_branching_chance() - curParams.min_branching_chance());
        if (dice(b_ch,1.0))
        {
            j.type = j.FORK;
            j.max_branching = 1 + (curParams.max_branching() - 1)*floor((float)rand()/RAND_MAX);
        }
        else j.type = j.MIDDLE;
    }
    voxels->set_occluder(j.pos,1.0);
    calc_light(j);    
    b->joints.push_back(j);
}
void TreeGenerator::new_branch(Branch *b, Joint &j, Segment &s, glm::vec3 &M, bool from_end)
{
    Branch *nb = curTree.branchHeaps[b->level + 1]->new_branch();
    nb->base_seg_n = from_end ? b->segments.size() : 0;
    nb->level = b->level + 1;
    nb->max_seg_count = curParams.max_segments();
    nb->base_r = curParams.base_r();
    Joint nj;
    nj.pos = j.pos;
    new_joint(nb, nj);
    glm::vec3 prev_dir = glm::normalize(s.end - s.begin);
    glm::vec3 rnd = rand_planar_dir();
    glm::vec3 N = glm::normalize(rnd - dot(rnd,prev_dir)*rnd); //Normal Vector
    glm::vec3 up = glm::vec3(0, 1, 0);

    glm::vec3 dir = curParams.dir_conserv() * prev_dir + curParams.spread() * N + curParams.phototrop() * M + curParams.gravitrop() * up;
    dir = from_end ? prev_dir : glm::normalize(dir);
    
    new_segment2(nb,dir,j.pos);
    calc_light(nb);
    calc_size(nb);
    nb->light += j.light;
    j.light = 0;

    j.childBranches.push_back(nb);
    test = nb;
}
void TreeGenerator::try_new_branch(Branch *b,Joint &j, Segment &s, bool from_end)
{   
    curParams.set_state(b->level+1);
    float feed = j.light;
    int bs = j.childBranches.size();
    if ((b->level<curParams.max_depth() - 1) 
       && ((j.type == j.FORK && &j != &(b->joints.back()) && (bs<j.max_branching)) || from_end)
       && dice(feed,exp(bs)*curParams.base_branch_feed()))
    {
        float occ = 0.0;
        glm::vec3 M = voxels->get_dir_to_bright_place(j.pos,&occ);
        if (dice(1, occ*curParams.branch_grow_decrease_q()))
        {
            new_branch(b,j,s,M,from_end);
        }
    }
    curParams.set_state(b->level);
}
inline float lfun(float v1,float v2,float a,float b)
{
    return (v1/v2)*(b-a) + a;
}
void TreeGenerator::set_seg_r(Branch *base, Segment &s)
{
    float b_st = curParams.base_r();
    if (base->level + 1 < curParams.max_depth())
        curParams.set_state(base->level+1);
    float b_end = curParams.base_r();
    curParams.set_state(base->level);
    int n = base->segments.size();
    s.rel_r_end = lfun(curParams.max_segments() - n-1,curParams.max_segments(),b_end,b_st);
    s.rel_r_begin = lfun(curParams.max_segments() - n-1,curParams.max_segments(),b_end,b_st);
}
void TreeGenerator::new_segment2(Branch *base, glm::vec3 &dir, glm::vec3 &pos)
{
    Segment s;
    float len = curParams.seg_len_mult();
    s.begin = pos;
    s.end = s.begin + len * dir;
    set_seg_r(base,s);
    Joint j;
    j.pos = s.end;
    base->segments.push_back(s);
    new_joint(base, j);
}
void TreeGenerator::new_segment(Branch *base, glm::vec3 &M)
{
    Joint &last = base->joints.back();
    Segment &last_s = base->segments.back();

    float sp = curParams.seg_spread();
    float bn = base->level*curParams.seg_bend()*powf((float)base->segments.size()/base->max_seg_count,curParams.seg_bend_pow());
    glm::vec3 prev_dir = glm::normalize(last_s.end - last_s.begin);
    glm::vec3 rnd = rand_dir();
    glm::vec3 N = glm::normalize(glm::cross(prev_dir, rnd)); //Normal Vector
    glm::vec3 dir = glm::normalize(glm::mix(prev_dir, N, sp));
    glm::vec3 up = glm::vec3(0, 1, 0);

    dir = curParams.seg_dir_conserv() * prev_dir + curParams.seg_spread ()* N + curParams.seg_phototrop() * M 
          + curParams.seg_gravitrop() * up - bn*up;
    dir = glm::normalize(dir);

    new_segment2(base,dir,last.pos);
}
void TreeGenerator::try_new_segment(Branch *base)
{
    float feed = base->joints.back().light;
    if ((base->segments.size()<base->max_seg_count) && dice(feed,curParams.base_seg_feed()))
    {
        float occ = 0.0;
        sum_feed[base->level] += feed;
        count_feed[base->level] += 1;
        glm::vec3 M = voxels->get_dir_to_bright_place(base->segments.back().end, &occ);
        if (dice(1,occ*curParams.segment_grow_decrease_q()))
        {
            new_segment(base,M);
        }
    }
}
bool comp(float *a, float *b)
{
    return *a>*b;
}
void TreeGenerator::distribute_feed(Branch *b)
{
    std::vector<float *> feeds;
    float *top = &(b->joints.back().light);
    //feeds.push_back(top);
    float min_light = 1.05;
    for (auto &j : b->joints)
    {
        if (j.childBranches.empty() && j.light > min_light)
            feeds.push_back(&(j.light));
        for (auto ch_b : j.childBranches)
        {
            if (ch_b->light > min_light)
                feeds.push_back( &(ch_b->light));
        }
    }
    if (false && b == test)
    {
        //fprintf(stderr,"id = %d)test feeds ",curTree.id);
        for (auto f:feeds)
        {
            //fprintf(stderr,"%f ",*f);
        }
        //fprintf(stderr,"\n");
    }
    std::sort(feeds.begin(),feeds.end(),comp);
    if (false && b == test)
    {
        //fprintf(stderr,"test feeds ");
        for (auto f:feeds)
        {
            //fprintf(stderr,"%f ",*f);
        }
        //fprintf(stderr,"\n");
    }
    float w = 1.0;
    float w_sum = 0.0;
    float w_decay = 1 - curParams.feed_distribution_d_weight();
    for (auto f: feeds)
    {
        (*f) = 1/w;
        w_sum += 1/w;
        w++;
        if (1/w<curParams.feed_distribution_min_weight())
            w = 1/curParams.feed_distribution_min_weight();
    }
    float top_bonus = is_branch_productive(b) ? w_sum*curParams.top_growth_bonus() / (1 - curParams.top_growth_bonus()) : 0;
    w_sum += top_bonus;
    *top += top_bonus;
    float mult = b->light/w_sum;
    if (false && b == test)
    {
        //fprintf(stderr,"results feeds ");
        for (auto f:feeds)
        {
            //fprintf(stderr,"%f ",*f);
        }
        //fprintf(stderr,"\n");
    }
    for (auto f: feeds)
    {
        (*f) *=mult;
    }
    if (false && b == test)
    {
        //fprintf(stderr,"mult feeds ");
        for (auto f:feeds)
        {
            //fprintf(stderr,"%uud %f ",f,*f);
        }
        //fprintf(stderr,"\n");
    }

}
void TreeGenerator::remove_branch(Branch *b)
{
    b->dead = true;
    bool save_trunk = false;
    Joint j1,j2;
    Segment s;
    if (b->joints.size()>=2 && b->segments.size() >=1)
    {
        save_trunk = true;
        auto jj = b->joints.begin();
        j1 = *jj;
        j2 = *(jj++);
        s = b->segments.front();
    }
    for (auto j : b->joints)
    {
        if (j.leaf)
            j.leaf->dead = true;
        for (auto br : j.childBranches)
            remove_branch(br);
    }
    b->joints.clear();
    b->segments.clear();
    if (save_trunk)
    {
        b->joints.push_back(j1);
        b->joints.push_back(j2);
        b->segments.push_back(s);
    }
}
void TreeGenerator::grow_branch(Branch *b, float feed)
{
    if (b->dead)
        return;
    curParams.set_state(b->level);
    float average_feed = b->light/(b->size+0.001);
    if (b->size && (b->base_seg_n == 0) && !dice(average_feed,curParams.branch_removal()))
    {
        remove_branch(b);
        //fprintf(stderr,"branch removed %f %f av_feed = %f\n",b->light, b->size, average_feed);
        return;
    }
    distribute_feed(b);

    auto j = b->joints.begin();
    j++;
    auto s = b->segments.begin();
    while ((j!=b->joints.end())&&(s!=b->segments.end()))
    {   
        bool from_end = (b->joints.size() == b->max_seg_count + 1) && (std::next(j) == b->joints.end()) && (j->childBranches.size() == 0);
        try_new_branch(b,*j,*s,from_end);
        j++;
        s++;
    }
    try_new_segment(b);
    recalculate_thickness(b);
    for (auto &j : b->joints)
    {
        for (auto br : j.childBranches)
        {
            grow_branch(br, 1);
        }
    }
}
float TreeGenerator::calc_light(Joint& j)
{
    float l = j.childBranches.empty() ? voxels->get_occlusion(j.pos) : 0;
    l = 50 - l;
    if (l<1)
        l = 1;
    j.light = l;
    for (auto b : j.childBranches)
    {
        l +=calc_light(b);
    }
    return l;
}
float TreeGenerator::calc_size(Joint& j)
{
    float sz = 0;
    for (auto b : j.childBranches)
    {
        sz += calc_size(b);
    }
    return sz;
}
float TreeGenerator::calc_size(Branch *b)
{
    float sz = b->segments.size()*powf(2,curParams.max_depth() - b->level);
    for (auto &j : b->joints)
    {
        sz +=calc_size(j);
    }
    b->size = sz;
    return sz;
}
float TreeGenerator::calc_light(Branch *b)
{
    if (b->dead)
        return 0;
    float l = 0.0;
    for (auto &j : b->joints)
    {
        l +=calc_light(j);
    }
    b->light = l;
    return l;
}
void TreeGenerator::recalculate_thickness(Branch *b)
{
    float *weights = new float[b->joints.size()+1];
    int i=b->joints.size()-1;
    weights[i+1] = 0.5;
    auto rev_it = b->joints.rbegin();
    while (rev_it != b->joints.rend())
    {
        weights[i] = weights[i+1] + rev_it->light;
        for (auto br: rev_it->childBranches)
        {
            weights[i]+=br->light;
        }
        i--;
        rev_it++;
    }
    i = 1;
    auto s_it = b->segments.begin();
    auto s_prev = b->segments.begin();
    auto j_it = b->joints.begin();
    j_it++;
    s_it++;
    s_prev->rel_r_begin = b->base_r;
    while (s_it != b->segments.end())
    {
        float sum_r = s_prev->rel_r_begin;
        s_it->rel_r_begin = sum_r*powf((weights[i]+ j_it->light)/weights[i-1],1/curParams.r_split_save_pow());
        s_prev->rel_r_end = s_it->rel_r_begin;
        for (auto br: j_it->childBranches)
            br->base_r = sum_r*powf(br->light/weights[i-1],1/curParams.r_split_save_pow());
        s_prev++;
        s_it++;
        i++;
        j_it++;
    }
    s_prev->rel_r_end = s_prev->rel_r_begin*powf(weights[i]/weights[i-1],1/curParams.r_split_save_pow());
    if (!b->joints.back().childBranches.empty())
    {
        b->joints.back().childBranches.front()->base_r = s_prev->rel_r_end;
    }
    delete[](weights);
}
LightVoxelsCube *TreeGenerator::create_light_voxels_cube(TreeStructureParameters params, glm::vec3 pos)
{
    float max_size = 0.0;
    float wm = 0.7;
    float hm = 0.55;
    for (int i=0;i<params.max_depth();i++)
    {
        params.set_state(i);
        max_size += params.max_segments()*params.seg_len_mult();
    }
    params.set_state(params.max_depth() - 1);
    return new LightVoxelsCube(pos + glm::vec3(0,0.5*max_size,0),glm::vec3(wm*max_size,hm*max_size,wm*max_size),
                               params.seg_len_mult(),params.light_precision());
}
void TreeGenerator::plant_tree(Tree &t, TreeStructureParameters params)
{
    for (int i=0;i<10;i++)
    {
        sum_feed[i] = 0;
        count_feed[i] = 0;
    }
    t.voxels = create_light_voxels_cube(params,t.pos);
    voxels = t.voxels;
    glm::vec3 sun_dir(-1,-1,-1);
    t.params = params;
    //voxels->set_directed_light(sun_dir,1.0);
    t.leaves = new LeafHeap();
    for (int i=0;i<params.max_depth();i++)
        t.branchHeaps.push_back(new BranchHeap());
    curParams.set_state(0);
    root = t.branchHeaps[0]->new_branch();
    root->level = 0;
    root->max_seg_count = curParams.max_segments();
    root->base_seg_n = 0;
    root->base_r = curParams.base_r();
    Joint j1,j2;
    Segment ts;
    j1.pos = t.pos;
    new_joint(root,j1);
    ts.begin = j1.pos;
    ts.end = j1.pos + glm::vec3(0,curParams.base_r(),0);
    j2.pos = ts.end;
    ts.rel_r_begin = 1;
    ts.rel_r_end = 1;
    
    new_joint(root,j2);
    root->segments.push_back(ts);
    curTree = t;
    curParams = params;
    t.root = root;
}
void TreeGenerator::grow_tree(Tree &t)
{
    curParams = t.params;
    curTree = t;
    root = t.root;
    if (t.branchHeaps.size()>1 && !t.branchHeaps[1]->branches.empty())
        test = &(t.branchHeaps[1]->branches.front());
    test = nullptr;
    voxels = t.voxels;
    for (int i=0;i<10;i++)
    {
        sum_feed[i] = 0;
        count_feed[i] = 0;
    }
    if (t.iter<curParams.growth_iterations() && root)
    {
        float feed = 1;
        seg_count = 0;
        feed = calc_light(root);
        root->light += 1000;
        calc_size(root);
        root->size = 0.01;
        grow_branch(root,feed);
        t.iter++;

        if (!(t.iter % 10))
        {   
            fprintf(stderr,"tree is growing  %d/%d iteration %d leaves\n",t.iter,curParams.growth_iterations(),t.leaves->leaves.size());
            fprintf(stderr,"sum feed %f\n",feed);
            fprintf(stderr,"average feed distribution:");
            for (int j=0;j<10;j++)
            {
                if (count_feed[j]<0.1)
                    fprintf(stderr," nan");
                else
                {
                    fprintf(stderr," %f",sum_feed[j]/count_feed[j]);
                }
                
            }
            fprintf(stderr,"\n");
            voxels->print_average_occlusion();
        }
    }
}
void TreeGenerator::leaf_to_model(Leaf &l, Model *m)
{
  if (l.dead)
    return;
  if (l.edges.size()<4)
  {
    //fprintf(stderr,"empty leaf\n");
    return;
  }
  glm::vec3 a = l.edges[0];
  glm::vec3 b = l.edges[1];
  glm::vec3 c = l.edges[2];
  glm::vec3 n = glm::normalize(glm::cross(a-b,c-b));
  std::vector<float> tex_c{1,0,0,0,0,1,1,1};
  int _b = m->positions.size()/3;
  for (int i=0;i<4;i++)
  {
    glm::vec3 v = l.edges[i];
    m->positions.push_back(v.x);
    m->positions.push_back(v.y);
    m->positions.push_back(v.z);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    m->colors.push_back(tex_c[2*i]);
    m->colors.push_back(tex_c[2*i+1]);
    m->colors.push_back(0);
    m->colors.push_back(1);
  }

  m->indices.push_back(_b);
  m->indices.push_back(_b+1);
  m->indices.push_back(_b+2);
  m->indices.push_back(_b+2);
  m->indices.push_back(_b+3);
  m->indices.push_back(_b);
}
void TreeGenerator::get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos)
{
  sv.ringsize = ring_size;
  glm::vec3 start = s.begin;
  glm::vec3 end   = s.end;
  glm::vec3 dir = end - start;
  
  //Get Some Normal Vector
  glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
  glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

  //Add the Correct Number of Indices
  glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)(2*PI/ring_size), dir);

  for (int i = 0; i < ring_size; i++)
  {
    VertexData vd;
    vd.pos = start + s.rel_r_begin*glm::vec3(n.x,n.y,n.z);
    vd.normal = n;
    vd.tex_coord.x = ((float)i)/ring_size;
    vd.tex_coord.y = rel_ring_pos;
    sv.bigRing.push_back(vd);
    n = r*n;
    ////fprintf(stderr,"size = %f\n",sv.bigRing.size());
  }
  ////fprintf(stderr,"size = %d %d\n",ringsize, sv.bigRing.size());
}
void TreeGenerator::get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos)
{
  sv.ringsize = ring_size;
  glm::vec3 start = s.begin;
  glm::vec3 end   = s.end;
  glm::vec3 dir = end - start;
  
  //Get Some Normal Vector
  glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
  glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

  //Add the Correct Number of Indices
  glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)(2*PI/ring_size), dir);

  for (int i = 0; i < ring_size; i++)
  {
    VertexData vd;
    vd.pos = end + s.rel_r_end*glm::vec3(n.x,n.y,n.z);
    vd.normal = n;
    vd.tex_coord.x = ((float)i)/ring_size;
    vd.tex_coord.y = rel_ring_pos;
    sv.smallRing.push_back(vd);
    n = r*n;
    ////fprintf(stderr,"size = %f\n",sv.bigRing.size());
  }
  ////fprintf(stderr,"size = %d %d\n",ringsize, sv.bigRing.size());
}
void pv(glm::vec3 pos)
{
  //fprintf(stderr,"(%f %f %f)",pos.x,pos.y,pos.z);
}
void pt(Model *m,int v)
{
  int t = m->indices[v];
  if (3*t + 3>m->positions.size())
  {
    //fprintf(stderr,"wrong index ind[%d] = %d",v,t);
  }
  else
  {
    glm::vec3 pos = glm::vec3(m->positions[3*t],m->positions[3*t+1],m->positions[3*t+2]);
    pv(pos);
  }
}
void TreeGenerator::seg_vertexes_to_model(SegmentVertexes &sv, Model *m)
{
  Model *h = m;
  int _b = h->positions.size()/3;
  if (sv.smallRing.size()<sv.ringsize || sv.bigRing.size()<sv.ringsize)
    return;
  ////fprintf(stderr,"_b = %d",_b);
  for (auto pos: sv.smallRing)
      {
        ////fprintf(stderr,"yy pos = %f %f %f\n",pos.pos.x,pos.pos.y,pos.pos.z);
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        h->colors.push_back(pos.tex_coord.x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(0.0);
        h->colors.push_back(1.0);
      }
      for (auto pos: sv.bigRing)
      {
        ////fprintf(stderr,"dd pos = %f %f %f\n",pos.pos.x,pos.pos.y,pos.pos.z);
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        h->colors.push_back(pos.tex_coord.x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(0.0);
        h->colors.push_back(1.0);
      }
  int ringsize = sv.ringsize;
  for(int i = 0; i < ringsize; i++)
  {
    //Bottom Triangle
    h->indices.push_back(_b+i);
    h->indices.push_back(_b+(i + 1)%(ringsize));
    h->indices.push_back(_b+i+ringsize);
    //Upper Triangle
    h->indices.push_back(_b+i+ringsize);
    h->indices.push_back(_b+(i + 1)%(ringsize));
    h->indices.push_back(_b+(i + 1)%(ringsize) + ringsize);
  }
}
void TreeGenerator::joint_to_model(Joint &j, Model *m, bool leaves)
{

}
void TreeGenerator::recursive_branch_to_model(Branch &b, Model *m, bool leaves)
{
    if (b.level<0)
        return;
    if (!leaves)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, curParams.max_depth() - b.level);
        for (auto &segment : b.segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b.segments.back(),vets.back(),ringsize,(float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m);
        }
    }
    else
    {
        for (auto &joint : b.joints)
        {
            if (joint.leaf)
                leaf_to_model(*(joint.leaf), m);
        }
    }
    for (auto &joint : b.joints)
    {
        for (auto branch : joint.childBranches)
            recursive_branch_to_model(*branch, m, leaves);
    }
}
void TreeGenerator::branch_to_model(Branch &b, Model *m, bool leaves)
{
    if (!leaves)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, curParams.max_depth() - b.level);
        for (auto &segment : b.segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b.segments.back(), vets.back(), ringsize, (float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m);
        }
    }
  else
  {
    for (auto &joint : b.joints)
        joint_to_model(joint,m,leaves);
  }
}
bool TreeGenerator::tree_to_model(Tree &t, bool leaves, bool debug)
{
    if (false && debug)
    {
        Model *m = new Model();
        for (int i = 0; i < 2/*t.branchHeaps.size()*/; i++)
        {
            for (auto &branch : t.branchHeaps[i]->branches)
            {
                if (!branch.dead)
                    branch_to_model(branch, m, leaves);
            }
        }  
        m->scale = glm::vec3(t.params.scale());
        t.models.push_back(m);
        BillboardCloud *cloud = new BillboardCloud(1, 1);
        t.billboardClouds.push_back(cloud);
        return true;
    }
    int sz = 1024;
    for (int level = 0; level < 3; level++)
    {
        Model *m = new Model();
        for (int i = 0; i < level; i++)
        {
            for (auto &branch : t.branchHeaps[i]->branches)
            {
                branch_to_model(branch, m, leaves);
            }
        }
        if (leaves)
        {
            for (auto &leaf : t.leaves->leaves)
            {
                leaf_to_model(leaf, m);
            }
        }
        m->scale = glm::vec3(t.params.scale());
        BillboardCloud *cloud = new BillboardCloud(sz, sz);
        cloud->prepare(t, level);
        t.billboardClouds.push_back(cloud);
        sz = sz * 2;
        t.models.push_back(m);
    }
    return true;
}
void TreeGenerator::segment_to_model(Segment &s, Model *m, bool leaves)
{
    
    glm::vec3 start = s.begin;
    glm::vec3 end   = s.end;
    glm::vec3 dir = end - start;
    int ringsize = 8;
    Model *h = m;
    //Get Some Normal Vector
    glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

    //Add the Correct Number of Indices
    glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)3.141/ringsize, dir);

    //Index Buffer
    int _b = h->positions.size()/3;

    //GL TRIANGLES
    for(int i = 0; i < ringsize; i++){
      //Bottom Triangle
      h->indices.push_back(_b+i*2+0);
      h->indices.push_back(_b+(i*2+2)%(2*ringsize));
      h->indices.push_back(_b+i*2+1);
      //Upper Triangle
      h->indices.push_back(_b+(i*2+2)%(2*ringsize));
      h->indices.push_back(_b+(i*2+3)%(2*ringsize));
      h->indices.push_back(_b+i*2+1);
    }

    for(int i = 0; i < ringsize; i++){

      h->positions.push_back(start.x + s.rel_r_begin*n.x);
      h->positions.push_back(start.y + s.rel_r_begin*n.y);
      h->positions.push_back(start.z + s.rel_r_begin*n.z);
      h->normals.push_back(n.x);
      h->normals.push_back(n.y);
      h->normals.push_back(n.z);
      n = r*n;

      h->positions.push_back(end.x + s.rel_r_end*n.x);
      h->positions.push_back(end.y + s.rel_r_end*n.y);
      h->positions.push_back(end.z + s.rel_r_end*n.z);
      h->normals.push_back(n.x);
      h->normals.push_back(n.y);
      h->normals.push_back(n.z);
      n = r*n;

    }
}
void LeafHeap::clear_removed()
{
}
void BranchHeap::clear_removed()
{

}
bool TreeGenerator::is_branch_productive(Branch *b)
{
    return (b->level < curParams.max_depth() - 1);
}
void TreeGenerator::create_tree(Tree &t, TreeStructureParameters params)
{
    plant_tree(t,params);
    while (t.iter < params.growth_iterations())
    {
        grow_tree(t);
    }
    Clusterizer cl;
    cl.set_branches(t,2);
    tree_to_model(t,false,false);
}