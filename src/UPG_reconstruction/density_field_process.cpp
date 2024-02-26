#include "density_field_process.h"

float* 
df::create_density_field(const std::vector<float>& model)
{
    std::srand(std::time(nullptr));

    struct df::VoxelGrid grid;
    grid.voxels = new df::Voxel [grid.dimension * grid.dimension * grid.dimension];

    float *density_field = new float [grid.dimension * grid.dimension * grid.dimension];
    
    glm::vec3 gridSize = grid.bounds[1] - grid.bounds[0];
    size_t numVoxels = grid.dimension * grid.dimension * grid.dimension;
    glm::vec3 voxelSize = gridSize / (float)grid.dimension;

    // std::cout << voxelSize.x << std::endl;

    for (int xi = 0; xi < grid.dimension; xi++)
    {
        for (int yi = 0; yi < grid.dimension; yi++)
        {
            for (int zi = 0; zi < grid.dimension; zi++)
            {
                glm::vec3 voxelCoord = glm::vec3(xi, yi, zi) * voxelSize + grid.bounds[0];
                df::Voxel Voxel = grid.voxels[(zi * grid.dimension + yi) * grid.dimension + xi];

                //  sampling every voxel and testing its points if they are inside a mesh
                
                for (int point_ind = 0; point_ind < samplePointsCount; point_ind++)
                {
                    Voxel.samplingPoints[point_ind] = glm::vec3(0 + (voxelSize.x - 0) * rand() / (float)RAND_MAX,
                                                                0 + (voxelSize.y - 0) * rand() / (float)RAND_MAX,
                                                                0 + (voxelSize.z - 0) * rand() / (float)RAND_MAX);
                    
                    //  check if point is inside
                    bool is_inside = true;

                    for (int i = 0; is_inside && i < model.size(); i += 3*FLOAT_PER_VERTEX)
                    {
                        //  Calculate angle between triangle norm and vector(sample point - triangle center) 

                        glm::vec3 norm = glm::vec3(model[i + 0 * FLOAT_PER_VERTEX + 3],
                                                    model[i + 0 * FLOAT_PER_VERTEX + 3 + 1],
                                                    model[i + 0 * FLOAT_PER_VERTEX + 3 + 2]);
                        
                        glm::vec3 tr_center = (glm::vec3(model[i+0*FLOAT_PER_VERTEX], model[i+0*FLOAT_PER_VERTEX+1], model[i+0*FLOAT_PER_VERTEX+2]) + 
                                                glm::vec3(model[i+1*FLOAT_PER_VERTEX], model[i+1*FLOAT_PER_VERTEX+1], model[i+1*FLOAT_PER_VERTEX+2]) + 
                                                glm::vec3(model[i+2*FLOAT_PER_VERTEX], model[i+2*FLOAT_PER_VERTEX+1], model[i+2*FLOAT_PER_VERTEX+2])) / 3.f;
                        
                        glm::vec3 p_vec = glm::normalize(voxelCoord + Voxel.samplingPoints[point_ind]);

                        float angle = glm::dot(norm, p_vec);
                        
                        is_inside = !(angle > 0);
                    }

                    if (is_inside)
                    {
                        // std::cout << 1;
                        Voxel.inside_count++;
                    }
                }

                //  Calc density in each voxel
                density_field[(zi * grid.dimension + yi) * grid.dimension + xi] = Voxel.inside_count / samplePointsCount;

                // std::cout << Voxel.inside_count << " " << samplePointsCount - Voxel.inside_count << std::endl;
            }
        }
    }

    return density_field;
}

void
df::erase(float* density, float trashold)
{
    for (size_t i = 0; i < density_dim * density_dim * density_dim; i++)
    {
        if (density[i] < trashold)
        {
            density[i] = 0.f;
        }
    }
}

void 
df::create_sdf(float* density)
{

}

float*
df::pipeline(const std::vector<float>& model)
{
    //  First step
    float* density = df::create_density_field(model);

    //  Second
    df::erase(density, 0.5);

    //  Third
    df::create_sdf(density);

    return nullptr;
}