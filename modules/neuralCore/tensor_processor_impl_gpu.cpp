#include <vector>
#include <memory>
#include <limits>
#include <cassert>
#include <chrono>
#include <array>

#include "vk_copy.h"
#include "vk_context.h"
#include "vk_images.h"

#include "tensor_processor_impl_gpu.h"
#include "include/TensorProcessorImpl_gpu_ubo.h"



std::shared_ptr<TensorProcessorImpl> CreateTensorProcessorImpl_GPU(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated) 
{ 
  auto pObj = std::make_shared<TensorProcessorImpl_GPU>(); 
  pObj->SetVulkanContext(a_ctx);
  pObj->InitVulkanObjects(a_ctx.device, a_ctx.physicalDevice, a_maxThreadsGenerated); 
  return pObj;
}

static uint32_t ComputeReductionSteps(uint32_t whole_size, uint32_t wg_size)
{
  uint32_t steps = 0;
  while (whole_size > 1)
  {
    steps++;
    whole_size = (whole_size + wg_size - 1) / wg_size;
  }
  return steps;
}

constexpr uint32_t KGEN_REDUCTION_LAST_STEP    = 16;

void TensorProcessorImpl_GPU::InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount) 
{
  physicalDevice = a_physicalDevice;
  device         = a_device;
  m_allCreatedPipelineLayouts.reserve(256);
  m_allCreatedPipelines.reserve(256);
  InitHelpers();
  InitBuffers(a_maxThreadsCount, true);
  InitKernels(".spv");
  AllocateAllDescriptorSets();

}

void TensorProcessorImpl_GPU::UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  const size_t maxAllowedSize = std::numeric_limits<uint32_t>::max();
  m_uboData.memory_size     = uint32_t( memory.size() );     assert( memory.size() < maxAllowedSize );
  m_uboData.memory_capacity = uint32_t( memory.capacity() ); assert( memory.capacity() < maxAllowedSize );
  a_pCopyEngine->UpdateBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
}

void TensorProcessorImpl_GPU::ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  a_pCopyEngine->ReadBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
  memory.resize(m_uboData.memory_size);
}

void TensorProcessorImpl_GPU::UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  if(memory.size() > 0)
    a_pCopyEngine->UpdateBuffer(m_vdata.memoryBuffer, 0, memory.data(), memory.size()*sizeof(float) );
}

void TensorProcessorImpl_GPU::UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{ 
}

void TensorProcessorImpl_GPU::cosCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, cosLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, cosPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::sinCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, sinLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sinPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::urandCmd(float *data, unsigned steps, Variable A, unsigned seed)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    unsigned int m_seed; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_seed = seed; 

  vkCmdPushConstants(m_currCmdBuffer, urandLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, urandPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::set_inputCmd(const float* data_in, unsigned offset, unsigned size)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_offset; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(size);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = size;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_offset = offset; 

  vkCmdPushConstants(m_currCmdBuffer, set_inputLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, set_inputPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::equalCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                         Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, equalLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, equalPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::minCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, minLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, minPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::logCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, logLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::logical_orCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                              Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, logical_orLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logical_orPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::not_equalCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                             Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, not_equalLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, not_equalPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::sqrtCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, sqrtLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sqrtPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::maxCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, maxLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, maxPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::powCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, powLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, powPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::greaterCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                           Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, greaterLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, greaterPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::dilateCmd(float *data, unsigned steps, unsigned x_size, unsigned x_dilate, unsigned y_size, unsigned y_dilate,
                                          unsigned z_size, unsigned z_dilate, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_x_size; 
    unsigned int m_x_dilate; 
    unsigned int m_y_size; 
    unsigned int m_y_dilate; 
    unsigned int m_z_size; 
    unsigned int m_z_dilate; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_x_size = x_size; 
  pcData.m_x_dilate = x_dilate; 
  pcData.m_y_size = y_size; 
  pcData.m_y_dilate = y_dilate; 
  pcData.m_z_size = z_size; 
  pcData.m_z_dilate = z_dilate; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, dilateLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, dilatePipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::outer_productCmd(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
                                      Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(A_len);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(B_len);
  
  pcData.m_sizeX  = A_len;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = B_len;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, outer_productLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, outer_productPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::padCmd(float *data, unsigned steps, unsigned step_size, unsigned left_pad, unsigned right_pad, 
                                       Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    unsigned int m_left_pad; 
    unsigned int m_right_pad; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m_left_pad = left_pad; 
  pcData.m_right_pad = right_pad; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, padLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, padPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::lessCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                        Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, lessLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, lessPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::max_pool_3D_diffCmd(float *data, int steps, int x_steps, int y_steps, int z_steps, 
                                                    int window_x, int window_y, int window_z, 
                                                    Variable A, Variable dLoss_dOutput, Variable dLoss_dInput)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_x_steps; 
    int m_window_x; 
    int m_window_y; 
    int m_window_z; 
    Variable m_A; 
    Variable m_dLoss_dOutput; 
    Variable m_dLoss_dInput; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(z_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(y_steps);
  
  pcData.m_sizeX  = z_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = y_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_x_steps = x_steps; 
  pcData.m_window_x = window_x; 
  pcData.m_window_y = window_y; 
  pcData.m_window_z = window_z; 
  pcData.m_A = A; 
  pcData.m_dLoss_dOutput = dLoss_dOutput; 
  pcData.m_dLoss_dInput = dLoss_dInput; 

  vkCmdPushConstants(m_currCmdBuffer, max_pool_3D_diffLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_3D_diffPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::less_equalCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                              Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, less_equalLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, less_equalPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::sumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, sumLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sumPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::copyCmd(float *data, unsigned steps, unsigned from, unsigned to, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_from; 
    unsigned int m_to; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_from = from; 
  pcData.m_to = to; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, copyLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::osumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, osumLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, osumPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::notCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, notLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, notPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::expCmd(float *data, unsigned steps, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, expLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, expPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::fillCmd(float *data, unsigned steps, Variable A, float val)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    Variable m_A; 
    float m_val; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A = A; 
  pcData.m_val = val; 

  vkCmdPushConstants(m_currCmdBuffer, fillLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, fillPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::get_outputCmd(float* data_out, unsigned offset, unsigned size)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_offset; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(size);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = size;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_offset = offset; 

  vkCmdPushConstants(m_currCmdBuffer, get_outputLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, get_outputPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::logical_andCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                               Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, logical_andLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logical_andPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::flipCmd(float *data, unsigned steps, unsigned flip_size, unsigned group_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_flip_size; 
    unsigned int m_group_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_flip_size = flip_size; 
  pcData.m_group_size = group_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, flipLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, flipPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::subCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, subLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, subPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::addCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, addLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, addPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::mulCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, mulLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, mulPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::greater_equalCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                                 Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, greater_equalLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, greater_equalPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::divCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, divLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, divPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::minimumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, minimumLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, minimumPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::max_pool_3DCmd(float *data, int steps, int x_steps, int y_steps, int z_steps, 
                                               int window_x, int window_y, int window_z, Variable A, Variable res)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_x_steps; 
    int m_window_x; 
    int m_window_y; 
    int m_window_z; 
    Variable m_A; 
    Variable m_res; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(z_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(y_steps);
  
  pcData.m_sizeX  = z_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = y_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_x_steps = x_steps; 
  pcData.m_window_x = window_x; 
  pcData.m_window_y = window_y; 
  pcData.m_window_z = window_z; 
  pcData.m_A = A; 
  pcData.m_res = res; 

  vkCmdPushConstants(m_currCmdBuffer, max_pool_3DLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_3DPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::maximumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, maximumLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, maximumPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::transposeCmd(float *data, unsigned steps, unsigned row_len, unsigned col_len, unsigned group_size, Variable A, Variable B)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    unsigned int m_group_size; 
    Variable m_A; 
    Variable m_B; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(row_len);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(col_len);
  
  pcData.m_sizeX  = row_len;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = col_len;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_group_size = group_size; 
  pcData.m_A = A; 
  pcData.m_B = B; 

  vkCmdPushConstants(m_currCmdBuffer, transposeLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, transposePipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::matmul_transposedCmd(float *data, unsigned A_col_len, unsigned B_col_len, 
                                                     unsigned A_row_len, Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 32;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_A_row_len; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(B_col_len);
  uint32_t sizeY  = uint32_t(A_col_len);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = B_col_len;
  pcData.m_sizeY  = A_col_len;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_A_row_len = A_row_len; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, matmul_transposedLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, matmul_transposedPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::max_poolCmd(float *data, int steps, int x_steps, int y_steps, int window_x, int window_y, 
                                            Variable A, Variable res)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_window_x; 
    int m_window_y; 
    Variable m_A; 
    Variable m_res; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(y_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(x_steps);
  
  pcData.m_sizeX  = y_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = x_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_window_x = window_x; 
  pcData.m_window_y = window_y; 
  pcData.m_A = A; 
  pcData.m_res = res; 

  vkCmdPushConstants(m_currCmdBuffer, max_poolLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_poolPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::whereCmd(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                         Variable A, Variable B, Variable C)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_total_size; 
    unsigned int m_step_size; 
    unsigned int m_group_size; 
    unsigned int m_Ai_mul; 
    Variable m_A; 
    Variable m_B; 
    Variable m_C; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_total_size = total_size; 
  pcData.m_step_size = step_size; 
  pcData.m_group_size = group_size; 
  pcData.m_Ai_mul = Ai_mul; 
  pcData.m_A = A; 
  pcData.m_B = B; 
  pcData.m_C = C; 

  vkCmdPushConstants(m_currCmdBuffer, whereLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, wherePipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::smax_diffCmd(float *data, unsigned steps, unsigned step_size, 
                                             Variable _output, Variable dLoss_dOutput, Variable dLoss_dInput)
{
  uint32_t blockSizeX = 256;
  uint32_t blockSizeY = 1;
  uint32_t blockSizeZ = 1;

  struct KernelArgsPC
  {
    unsigned int m_step_size; 
    Variable m__output; 
    Variable m_dLoss_dOutput; 
    Variable m_dLoss_dInput; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(steps);
  uint32_t sizeY  = uint32_t(1);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = steps;
  pcData.m_sizeY  = 1;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_step_size = step_size; 
  pcData.m__output = _output; 
  pcData.m_dLoss_dOutput = dLoss_dOutput; 
  pcData.m_dLoss_dInput = dLoss_dInput; 

  vkCmdPushConstants(m_currCmdBuffer, smax_diffLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, smax_diffPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::conv3dCmd(float *data, int steps, int x_steps, int y_steps, int z_steps, int stride, 
                                          int in_channels, int out_channels, Variable A, Variable kernel, Variable res)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_x_steps; 
    int m_stride; 
    int m_in_channels; 
    int m_out_channels; 
    Variable m_A; 
    Variable m_kernel; 
    Variable m_res; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(z_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(y_steps);
  
  pcData.m_sizeX  = z_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = y_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_x_steps = x_steps; 
  pcData.m_stride = stride; 
  pcData.m_in_channels = in_channels; 
  pcData.m_out_channels = out_channels; 
  pcData.m_A = A; 
  pcData.m_kernel = kernel; 
  pcData.m_res = res; 

  vkCmdPushConstants(m_currCmdBuffer, conv3dLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, conv3dPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::conv2dCmd(float *data, int steps, int x_steps, int y_steps, int stride, int in_channels, 
                                          int out_channels, Variable A, Variable kernel, Variable res)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_stride; 
    int m_in_channels; 
    int m_out_channels; 
    Variable m_A; 
    Variable m_kernel; 
    Variable m_res; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(y_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(x_steps);
  
  pcData.m_sizeX  = y_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = x_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_stride = stride; 
  pcData.m_in_channels = in_channels; 
  pcData.m_out_channels = out_channels; 
  pcData.m_A = A; 
  pcData.m_kernel = kernel; 
  pcData.m_res = res; 

  vkCmdPushConstants(m_currCmdBuffer, conv2dLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, conv2dPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}

void TensorProcessorImpl_GPU::max_pool_diffCmd(float *data, int steps, int x_steps, int y_steps, int window_x, int window_y,
                                                 Variable A, Variable dLoss_dOutput, Variable dLoss_dInput)
{
  uint32_t blockSizeX = 8;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 8;

  struct KernelArgsPC
  {
    int m_window_x; 
    int m_window_y; 
    Variable m_A; 
    Variable m_dLoss_dOutput; 
    Variable m_dLoss_dInput; 
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(y_steps);
  uint32_t sizeY  = uint32_t(steps);
  uint32_t sizeZ  = uint32_t(x_steps);
  
  pcData.m_sizeX  = y_steps;
  pcData.m_sizeY  = steps;
  pcData.m_sizeZ  = x_steps;
  pcData.m_tFlags = m_currThreadFlags;
  pcData.m_window_x = window_x; 
  pcData.m_window_y = window_y; 
  pcData.m_A = A; 
  pcData.m_dLoss_dOutput = dLoss_dOutput; 
  pcData.m_dLoss_dInput = dLoss_dInput; 

  vkCmdPushConstants(m_currCmdBuffer, max_pool_diffLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_diffPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
 
}


void TensorProcessorImpl_GPU::copyKernelFloatCmd(uint32_t length)
{
  uint32_t blockSizeX = MEMCPY_BLOCK_SIZE;

  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyKernelFloatPipeline);
  vkCmdPushConstants(m_currCmdBuffer, copyKernelFloatLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &length);
  vkCmdDispatch(m_currCmdBuffer, (length + blockSizeX - 1) / blockSizeX, 1, 1);
}

VkBufferMemoryBarrier TensorProcessorImpl_GPU::BarrierForClearFlags(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_TRANSFER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

VkBufferMemoryBarrier TensorProcessorImpl_GPU::BarrierForSingleBuffer(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

void TensorProcessorImpl_GPU::BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum)
{
  for(uint32_t i=0; i<a_buffersNum;i++)
  {
    a_outBarriers[i].sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
    a_outBarriers[i].pNext               = NULL;
    a_outBarriers[i].srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
    a_outBarriers[i].dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
    a_outBarriers[i].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].buffer              = a_inBuffers[i];
    a_outBarriers[i].offset              = 0;
    a_outBarriers[i].size                = VK_WHOLE_SIZE;
  }
}

void TensorProcessorImpl_GPU::get_outputCmd(VkCommandBuffer a_commandBuffer, float* out, unsigned offset, unsigned size)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
    vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, get_outputLayout, 0, 1, &m_allGeneratedDS[0], 0, nullptr);
  get_outputCmd(out, offset, size);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);

}

void TensorProcessorImpl_GPU::set_inputCmd(VkCommandBuffer a_commandBuffer, const float* in, unsigned offset, unsigned size)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
    vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, set_inputLayout, 0, 1, &m_allGeneratedDS[1], 0, nullptr);
  set_inputCmd(in, offset, size);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);

}

void TensorProcessorImpl_GPU::processCmd(VkCommandBuffer a_commandBuffer, const nn::TensorProgram &program)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
    unsigned data_size = memory.size();
  unsigned constexpr maxWorkGroupSizeY = (1<<16) - 1;
  unsigned constexpr maxWorkGroupSizeZ = (1<<16) - 1;
  #if DEBUG
  {
    printf("data [");
    for (int i=0;i<data_size;i++)
      printf("%8d ", i);
    printf("]\n");
  }
  #endif

  for (int i = 0; i < program.commands.size(); i++)
  {
    auto t1 = std::chrono::high_resolution_clock::now();

    Variable A = program.vars[program.commands[i].args[0]];
    Variable B = program.vars[program.commands[i].args[1]];
    Variable C = program.vars[program.commands[i].args[2]];

    unsigned arg0 = program.commands[i].args[3];
    unsigned arg1 = program.commands[i].args[4];
    unsigned arg2 = program.commands[i].args[5];
    unsigned arg3 = program.commands[i].args[6];
    unsigned arg4 = program.commands[i].args[7];

    switch (program.commands[i].type)
    {
    case nn::TensorProgram::NOOP:
      break;
    case nn::TensorProgram::ADD:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, addLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  addCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::MUL:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, mulLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  mulCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::SUB:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, subLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  subCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::DIV:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, divLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  divCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::GREATER:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, greaterLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  greaterCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::LESS:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, lessLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  lessCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::EQUAL:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, equalLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  equalCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::GE:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, greater_equalLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  greater_equalCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::LE:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, less_equalLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  less_equalCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::NE:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, not_equalLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  not_equalCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::OR:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logical_orLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  logical_orCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::AND:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logical_andLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  logical_andCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::WHERE:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, whereLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  whereCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;      
    case nn::TensorProgram::MIN:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, minLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  minCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;     
    case nn::TensorProgram::MAX:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, maxLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  maxCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;     
    case nn::TensorProgram::POW:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, powLayout, 0, 1, &m_allGeneratedDS[2], 0, nullptr);
  powCmd(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;     
    case nn::TensorProgram::EXP:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, expLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  expCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::SQRT:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sqrtLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  sqrtCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::SIN:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sinLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  sinCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::COS:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, cosLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  cosCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::LOG:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, logLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  logCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::NOT:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, notLayout, 0, 1, &m_allGeneratedDS[3], 0, nullptr);
  notCmd(memory.data(), A.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::SUM:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, sumLayout, 0, 1, &m_allGeneratedDS[4], 0, nullptr);
  sumCmd(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::O_SUM:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, osumLayout, 0, 1, &m_allGeneratedDS[4], 0, nullptr);
  osumCmd(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::MINIMUM:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, minimumLayout, 0, 1, &m_allGeneratedDS[4], 0, nullptr);
  minimumCmd(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::MAXIMUM:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, maximumLayout, 0, 1, &m_allGeneratedDS[4], 0, nullptr);
  maximumCmd(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::MATMUL_T:
      #if DEBUG
      if (B.Dim == 2 && B.sizes[1] > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MATMUL_T workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", B.sizes[1]);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, matmul_transposedLayout, 0, 1, &m_allGeneratedDS[5], 0, nullptr);
  matmul_transposedCmd(memory.data(), A.sizes[1], B.Dim == 2 ? B.sizes[1] : 1, A.sizes[0], A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::MOV:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyLayout, 0, 1, &m_allGeneratedDS[6], 0, nullptr);
  copyCmd(memory.data(), A.total_size, 0, 0, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::FILL:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, fillLayout, 0, 1, &m_allGeneratedDS[7], 0, nullptr);
  fillCmd(memory.data(), C.total_size, C, *((float *)(&arg0)));
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::COPY:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyLayout, 0, 1, &m_allGeneratedDS[8], 0, nullptr);
  copyCmd(memory.data(), arg2, arg0, arg1, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::TRANSP:
    {
      unsigned transp_dim = arg0;
      unsigned group_size = 1;
      for (unsigned d=0;d<transp_dim;d++)
        group_size *= A.sizes[d];
      unsigned steps = A.total_size/(A.sizes[transp_dim]*A.sizes[transp_dim+1]*group_size);
      unsigned row_len = A.sizes[transp_dim];
      unsigned col_len = A.sizes[transp_dim+1];
      #if DEBUG
      if (row_len > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: TRANSP workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", row_len);
      if (col_len > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: TRANSP workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", col_len);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, transposeLayout, 0, 1, &m_allGeneratedDS[9], 0, nullptr);
  transposeCmd(memory.data(), steps, row_len, col_len, group_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::OUTER_P:
      #if DEBUG
      if (A.sizes[0] > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: OUTER_P workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]);
      if (B.sizes[0] > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: OUTER_P workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", B.sizes[0]);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, outer_productLayout, 0, 1, &m_allGeneratedDS[10], 0, nullptr);
  outer_productCmd(memory.data(), A.total_size/A.sizes[0], A.sizes[0], B.sizes[0], A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::SMAX_D:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, smax_diffLayout, 0, 1, &m_allGeneratedDS[11], 0, nullptr);
  smax_diffCmd(memory.data(), A.sizes[A.Dim-1], A.total_size/A.sizes[A.Dim-1], A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::PAD:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, padLayout, 0, 1, &m_allGeneratedDS[12], 0, nullptr);
  padCmd(memory.data(), arg0, arg1, arg2, arg3, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::CONV_2D:
    {
      int in_channels = B.sizes[2] > 0 ? B.sizes[2] : 1;
      int out_channels = B.sizes[3] > 0 ? B.sizes[3] : 1;
      int images = A.total_size/(A.sizes[0]*A.sizes[1]*in_channels);
      int steps = images * out_channels;
      int x_steps = C.sizes[0];
      int y_steps = C.sizes[1];
      int stride = arg0;
      #if DEBUG
      if (x_steps > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: CONV_2D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", x_steps);
      if (y_steps > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: CONV_2D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", y_steps);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, conv2dLayout, 0, 1, &m_allGeneratedDS[13], 0, nullptr);
  conv2dCmd(memory.data(), steps, x_steps, y_steps, stride, in_channels, out_channels, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::FLIP:
    {
      unsigned flip_dim = arg0;
      unsigned group_size = 1;
      for (unsigned d=0;d<flip_dim;d++)
        group_size *= A.sizes[d];
      unsigned flip_size = A.sizes[flip_dim];
      unsigned steps = A.total_size/(flip_size*group_size);
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, flipLayout, 0, 1, &m_allGeneratedDS[14], 0, nullptr);
  flipCmd(memory.data(), steps, flip_size, group_size, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::MPOOL:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_poolLayout, 0, 1, &m_allGeneratedDS[15], 0, nullptr);
  max_poolCmd(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, window_x, window_y, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::MPOOL_D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_diffLayout, 0, 1, &m_allGeneratedDS[16], 0, nullptr);
  max_pool_diffCmd(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, window_x, window_y, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::DILATE:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, fillLayout, 0, 1, &m_allGeneratedDS[17], 0, nullptr);
  fillCmd(memory.data(), C.total_size, C, 0.0f);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, dilateLayout, 0, 1, &m_allGeneratedDS[18], 0, nullptr);
  dilateCmd(memory.data(), A.total_size, A.sizes[0], arg0, A.Dim>1?A.sizes[1]:1, arg1, A.Dim>2?A.sizes[2]:1, arg2, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::URAND:
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, urandLayout, 0, 1, &m_allGeneratedDS[19], 0, nullptr);
  urandCmd(memory.data(), C.total_size, C, arg0 == 0 ? rand() : arg0);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
      break;
    case nn::TensorProgram::CONV_3D:
    {
      int in_channels = B.sizes[3] > 0 ? B.sizes[3] : 1;
      int out_channels = B.sizes[4] > 0 ? B.sizes[4] : 1;
      int images = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]*in_channels);
      int steps = images * out_channels;
      int x_steps = C.sizes[0];
      int y_steps = C.sizes[1];
      int z_steps = C.sizes[2];
      int stride = arg0;
      #if DEBUG
      if (x_steps > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: CONV_3D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", x_steps);
      if (y_steps > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: CONV_3D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", y_steps);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, conv3dLayout, 0, 1, &m_allGeneratedDS[20], 0, nullptr);
  conv3dCmd(memory.data(), steps, x_steps, y_steps, z_steps, stride, in_channels, out_channels, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::MPOOL_3D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned window_z = arg2;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_3DLayout, 0, 1, &m_allGeneratedDS[21], 0, nullptr);
  max_pool_3DCmd(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, A.sizes[2]/window_z, window_x, window_y, window_z, A, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    case nn::TensorProgram::MPOOL_3D_D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned window_z = arg2;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, max_pool_3D_diffLayout, 0, 1, &m_allGeneratedDS[22], 0, nullptr);
  max_pool_3D_diffCmd(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, A.sizes[2]/window_z, window_x, window_y, window_z, A, B, C);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
    }
      break;
    default:
      break;
    }
    #if DEBUG
    {
      printf("data [");
      for (int i=0;i<data_size;i++)
        printf("%8.4f ", memory.data()[i]);
      printf("]\n");
    }
    #endif
    /*
    printf("cmd %d %s, C = [", i, nn::TensorProgram::cmd_properties[program.commands[i].type].name.c_str());
    bool has_nan = false;
    long double min = 1e15;
    long double max = -1e15;
    long double a = 0.0;
    long double aa = 0.0;
    for (int k=0;k<C.total_size;k++)
    {
      min = std::min(min, (long double)memory[C.offset+k]);
      max = std::max(max, (long double)memory[C.offset+k]);
      a += memory[C.offset+k];
      aa += std::abs(memory[C.offset+k]);
      if (memory[C.offset+k] != memory[C.offset+k])
        has_nan = true;
    }
    printf("min max av abs_av has_nan %f %f %f %f %d]\n",(float)min, (float)max,(float)a/C.total_size,(float)aa/C.total_size,(int)has_nan);*/

    auto t2 = std::chrono::high_resolution_clock::now();
    float total_time_us = 1e-3 * std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    _stat_time_cmd_id[program.commands[i].type] += total_time_us;
    _stat_time_cmd_num[i] += total_time_us;
  }
  _stat_execution_times++;

}



void TensorProcessorImpl_GPU::get_output(float* out, unsigned offset, unsigned size)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(get_output_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer outGPU = vk_utils::createBuffer(device, size*sizeof(float ), outFlags);
  buffers.push_back(outGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeget_output.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeget_output.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_get_output(outGPU, 0, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeget_output.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeget_output.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeget_output.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    get_outputCmd(commandBuffer, out, offset, size);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeget_output.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(outGPU, 0, out, size*sizeof(float ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeget_output.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, outGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeget_output.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void TensorProcessorImpl_GPU::set_input(const float* in, unsigned offset, unsigned size)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(set_input_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer inGPU = vk_utils::createBuffer(device, size*sizeof(const float ), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  buffers.push_back(inGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeset_input.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeset_input.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_set_input(inGPU, 0, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeset_input.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  pCopyHelper->UpdateBuffer(inGPU, 0, in, size*sizeof(const float )); 
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeset_input.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeset_input.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    set_inputCmd(commandBuffer, in, offset, size);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeset_input.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeset_input.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  //vkDestroyBuffer(device, inGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeset_input.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void TensorProcessorImpl_GPU::process(const nn::TensorProgram &program)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(process_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeprocess.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeprocess.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_process(); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeprocess.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeprocess.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeprocess.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    processCmd(commandBuffer, program);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeprocess.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeprocess.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeprocess.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void TensorProcessorImpl_GPU::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  vk_utils::ExecTime res = {};
  if(std::string(a_funcName) == "get_output" || std::string(a_funcName) == "get_outputBlock")
    res = m_exTimeget_output;
  if(std::string(a_funcName) == "set_input" || std::string(a_funcName) == "set_inputBlock")
    res = m_exTimeset_input;
  if(std::string(a_funcName) == "process" || std::string(a_funcName) == "processBlock")
    res = m_exTimeprocess;
  a_out[0] = res.msExecuteOnGPU;
  a_out[1] = res.msCopyToGPU;
  a_out[2] = res.msCopyFromGPU;
  a_out[3] = res.msAPIOverhead;             
}

