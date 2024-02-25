#include <vector>
#include <array>
#include <memory>
#include <limits>
#include <cassert>
#include "vk_copy.h"
#include "vk_context.h"

#include "tensor_processor_impl_gpu.h"
#include "include/TensorProcessorImpl_gpu_ubo.h"

static uint32_t ComputeReductionAuxBufferElements(uint32_t whole_size, uint32_t wg_size)
{
  uint32_t sizeTotal = 0;
  while (whole_size > 1)
  {
    whole_size  = (whole_size + wg_size - 1) / wg_size;
    sizeTotal  += std::max<uint32_t>(whole_size, 1);
  }
  return sizeTotal;
}

VkBufferUsageFlags TensorProcessorImpl_GPU::GetAdditionalFlagsForUBO() const
{
  return VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
}

uint32_t TensorProcessorImpl_GPU::GetDefaultMaxTextures() const { return 256; }

void TensorProcessorImpl_GPU::MakeComputePipelineAndLayout(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout* pPipelineLayout, VkPipeline* pPipeline)
{
  VkPipelineShaderStageCreateInfo shaderStageInfo = {};
  shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
  shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;

  auto shaderCode   = vk_utils::readSPVFile(a_shaderPath);
  auto shaderModule = vk_utils::createShaderModule(device, shaderCode);

  shaderStageInfo.module              = shaderModule;
  shaderStageInfo.pName               = a_mainName;
  shaderStageInfo.pSpecializationInfo = a_specInfo;

  VkPushConstantRange pcRange = {};
  pcRange.stageFlags = shaderStageInfo.stage;
  pcRange.offset     = 0;
  pcRange.size       = 128; // at least 128 bytes for push constants for all Vulkan implementations

  VkPipelineLayoutCreateInfo pipelineLayoutInfo = {};
  pipelineLayoutInfo.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
  pipelineLayoutInfo.pushConstantRangeCount = 1;
  pipelineLayoutInfo.pPushConstantRanges    = &pcRange;
  pipelineLayoutInfo.pSetLayouts            = &a_dsLayout;
  pipelineLayoutInfo.setLayoutCount         = 1;
   
  VkResult res = vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, pPipelineLayout);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreatePipelineLayout have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelineLayouts.push_back(*pPipelineLayout);

  VkComputePipelineCreateInfo pipelineInfo = {};
  pipelineInfo.sType              = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
  pipelineInfo.flags              = 0;
  pipelineInfo.stage              = shaderStageInfo;
  pipelineInfo.layout             = (*pPipelineLayout);
  pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;
  res = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, pPipeline);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreateComputePipelines have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelines.push_back(*pPipeline);

  if (shaderModule != VK_NULL_HANDLE)
    vkDestroyShaderModule(device, shaderModule, VK_NULL_HANDLE);
}

void TensorProcessorImpl_GPU::MakeComputePipelineOnly(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout pipelineLayout, VkPipeline* pPipeline)
{
  VkPipelineShaderStageCreateInfo shaderStageInfo = {};
  shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
  shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;

  auto shaderCode   = vk_utils::readSPVFile(a_shaderPath);
  auto shaderModule = vk_utils::createShaderModule(device, shaderCode);

  shaderStageInfo.module              = shaderModule;
  shaderStageInfo.pName               = a_mainName;
  shaderStageInfo.pSpecializationInfo = a_specInfo;

  VkComputePipelineCreateInfo pipelineInfo = {};
  pipelineInfo.sType              = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
  pipelineInfo.flags              = 0;
  pipelineInfo.stage              = shaderStageInfo;
  pipelineInfo.layout             = pipelineLayout;
  pipelineInfo.basePipelineHandle = VK_NULL_HANDLE;
  VkResult res = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, pPipeline);
  if(res != VK_SUCCESS)
  {
    std::string errMsg = vk_utils::errorString(res);
    std::cout << "[ShaderError]: vkCreateComputePipelines have failed for '" << a_shaderPath << "' with '" << errMsg.c_str() << "'" << std::endl;
  }
  else
    m_allCreatedPipelines.push_back(*pPipeline);

  if (shaderModule != VK_NULL_HANDLE)
    vkDestroyShaderModule(device, shaderModule, VK_NULL_HANDLE);
}


TensorProcessorImpl_GPU::~TensorProcessorImpl_GPU()
{
  for(size_t i=0;i<m_allCreatedPipelines.size();i++)
    vkDestroyPipeline(device, m_allCreatedPipelines[i], nullptr);
  for(size_t i=0;i<m_allCreatedPipelineLayouts.size();i++)
    vkDestroyPipelineLayout(device, m_allCreatedPipelineLayouts[i], nullptr);

  vkDestroyDescriptorSetLayout(device, cosDSLayout, nullptr);
  cosDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, sinDSLayout, nullptr);
  sinDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, urandDSLayout, nullptr);
  urandDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, set_inputDSLayout, nullptr);
  set_inputDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, equalDSLayout, nullptr);
  equalDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, minDSLayout, nullptr);
  minDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, logDSLayout, nullptr);
  logDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, logical_orDSLayout, nullptr);
  logical_orDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, not_equalDSLayout, nullptr);
  not_equalDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, sqrtDSLayout, nullptr);
  sqrtDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, maxDSLayout, nullptr);
  maxDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, powDSLayout, nullptr);
  powDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, greaterDSLayout, nullptr);
  greaterDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, dilateDSLayout, nullptr);
  dilateDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, outer_productDSLayout, nullptr);
  outer_productDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, padDSLayout, nullptr);
  padDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, lessDSLayout, nullptr);
  lessDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, max_pool_3D_diffDSLayout, nullptr);
  max_pool_3D_diffDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, less_equalDSLayout, nullptr);
  less_equalDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, sumDSLayout, nullptr);
  sumDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, copyDSLayout, nullptr);
  copyDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, osumDSLayout, nullptr);
  osumDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, notDSLayout, nullptr);
  notDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, expDSLayout, nullptr);
  expDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, fillDSLayout, nullptr);
  fillDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, get_outputDSLayout, nullptr);
  get_outputDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, logical_andDSLayout, nullptr);
  logical_andDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, flipDSLayout, nullptr);
  flipDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, subDSLayout, nullptr);
  subDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, addDSLayout, nullptr);
  addDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, mulDSLayout, nullptr);
  mulDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, greater_equalDSLayout, nullptr);
  greater_equalDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, divDSLayout, nullptr);
  divDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, minimumDSLayout, nullptr);
  minimumDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, max_pool_3DDSLayout, nullptr);
  max_pool_3DDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, maximumDSLayout, nullptr);
  maximumDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, transposeDSLayout, nullptr);
  transposeDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, matmul_transposedDSLayout, nullptr);
  matmul_transposedDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, max_poolDSLayout, nullptr);
  max_poolDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, whereDSLayout, nullptr);
  whereDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, smax_diffDSLayout, nullptr);
  smax_diffDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, conv3dDSLayout, nullptr);
  conv3dDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, conv2dDSLayout, nullptr);
  conv2dDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorSetLayout(device, max_pool_diffDSLayout, nullptr);
  max_pool_diffDSLayout = VK_NULL_HANDLE;
  vkDestroyDescriptorPool(device, m_dsPool, NULL); m_dsPool = VK_NULL_HANDLE;

 
  vkDestroyBuffer(device, m_classDataBuffer, nullptr);

  vkDestroyBuffer(device, m_vdata.memoryBuffer, nullptr);
  FreeAllAllocations(m_allMems);
}

void TensorProcessorImpl_GPU::InitHelpers()
{
  vkGetPhysicalDeviceProperties(physicalDevice, &m_devProps);
}


void TensorProcessorImpl_GPU::InitKernel_cos(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_cos.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  cosDSLayout = CreatecosDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, cosDSLayout, &cosLayout, &cosPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_sin(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_sin.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  sinDSLayout = CreatesinDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, sinDSLayout, &sinLayout, &sinPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_urand(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_urand.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  urandDSLayout = CreateurandDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, urandDSLayout, &urandLayout, &urandPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_set_input(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_set_input.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  set_inputDSLayout = Createset_inputDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, set_inputDSLayout, &set_inputLayout, &set_inputPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_equal(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_equal.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  equalDSLayout = CreateequalDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, equalDSLayout, &equalLayout, &equalPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_min(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_min.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  minDSLayout = CreateminDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, minDSLayout, &minLayout, &minPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_log(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_log.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  logDSLayout = CreatelogDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, logDSLayout, &logLayout, &logPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_logical_or(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_logical_or.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  logical_orDSLayout = Createlogical_orDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, logical_orDSLayout, &logical_orLayout, &logical_orPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_not_equal(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_not_equal.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  not_equalDSLayout = Createnot_equalDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, not_equalDSLayout, &not_equalLayout, &not_equalPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_sqrt(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_sqrt.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  sqrtDSLayout = CreatesqrtDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, sqrtDSLayout, &sqrtLayout, &sqrtPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_max(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_max.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  maxDSLayout = CreatemaxDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, maxDSLayout, &maxLayout, &maxPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_pow(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_pow.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  powDSLayout = CreatepowDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, powDSLayout, &powLayout, &powPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_greater(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_greater.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  greaterDSLayout = CreategreaterDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, greaterDSLayout, &greaterLayout, &greaterPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_dilate(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_dilate.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  dilateDSLayout = CreatedilateDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, dilateDSLayout, &dilateLayout, &dilatePipeline);
}

void TensorProcessorImpl_GPU::InitKernel_outer_product(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_outer_product.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  outer_productDSLayout = Createouter_productDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, outer_productDSLayout, &outer_productLayout, &outer_productPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_pad(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_pad.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  padDSLayout = CreatepadDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, padDSLayout, &padLayout, &padPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_less(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_less.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  lessDSLayout = CreatelessDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, lessDSLayout, &lessLayout, &lessPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_max_pool_3D_diff(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_max_pool_3D_diff.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  max_pool_3D_diffDSLayout = Createmax_pool_3D_diffDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, max_pool_3D_diffDSLayout, &max_pool_3D_diffLayout, &max_pool_3D_diffPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_less_equal(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_less_equal.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  less_equalDSLayout = Createless_equalDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, less_equalDSLayout, &less_equalLayout, &less_equalPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_sum(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_sum.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  sumDSLayout = CreatesumDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, sumDSLayout, &sumLayout, &sumPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_copy(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_copy.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  copyDSLayout = CreatecopyDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, copyDSLayout, &copyLayout, &copyPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_osum(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_osum.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  osumDSLayout = CreateosumDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, osumDSLayout, &osumLayout, &osumPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_not(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_not.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  notDSLayout = CreatenotDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, notDSLayout, &notLayout, &notPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_exp(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_exp.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  expDSLayout = CreateexpDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, expDSLayout, &expLayout, &expPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_fill(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_fill.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  fillDSLayout = CreatefillDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, fillDSLayout, &fillLayout, &fillPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_get_output(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_get_output.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  get_outputDSLayout = Createget_outputDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, get_outputDSLayout, &get_outputLayout, &get_outputPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_logical_and(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_logical_and.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  logical_andDSLayout = Createlogical_andDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, logical_andDSLayout, &logical_andLayout, &logical_andPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_flip(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_flip.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  flipDSLayout = CreateflipDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, flipDSLayout, &flipLayout, &flipPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_sub(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_sub.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  subDSLayout = CreatesubDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, subDSLayout, &subLayout, &subPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_add(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_add.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  addDSLayout = CreateaddDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, addDSLayout, &addLayout, &addPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_mul(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_mul.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  mulDSLayout = CreatemulDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, mulDSLayout, &mulLayout, &mulPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_greater_equal(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_greater_equal.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  greater_equalDSLayout = Creategreater_equalDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, greater_equalDSLayout, &greater_equalLayout, &greater_equalPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_div(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_div.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  divDSLayout = CreatedivDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, divDSLayout, &divLayout, &divPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_minimum(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_minimum.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  minimumDSLayout = CreateminimumDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, minimumDSLayout, &minimumLayout, &minimumPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_max_pool_3D(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_max_pool_3D.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  max_pool_3DDSLayout = Createmax_pool_3DDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, max_pool_3DDSLayout, &max_pool_3DLayout, &max_pool_3DPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_maximum(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_maximum.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  maximumDSLayout = CreatemaximumDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, maximumDSLayout, &maximumLayout, &maximumPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_transpose(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_transpose.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  transposeDSLayout = CreatetransposeDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, transposeDSLayout, &transposeLayout, &transposePipeline);
}

void TensorProcessorImpl_GPU::InitKernel_matmul_transposed(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel2D_matmul_transposed.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  matmul_transposedDSLayout = Creatematmul_transposedDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, matmul_transposedDSLayout, &matmul_transposedLayout, &matmul_transposedPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_max_pool(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_max_pool.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  max_poolDSLayout = Createmax_poolDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, max_poolDSLayout, &max_poolLayout, &max_poolPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_where(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_where.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  whereDSLayout = CreatewhereDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, whereDSLayout, &whereLayout, &wherePipeline);
}

void TensorProcessorImpl_GPU::InitKernel_smax_diff(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel1D_smax_diff.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  smax_diffDSLayout = Createsmax_diffDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, smax_diffDSLayout, &smax_diffLayout, &smax_diffPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_conv3d(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_conv3d.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  conv3dDSLayout = Createconv3dDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, conv3dDSLayout, &conv3dLayout, &conv3dPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_conv2d(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_conv2d.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  conv2dDSLayout = Createconv2dDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, conv2dDSLayout, &conv2dLayout, &conv2dPipeline);
}

void TensorProcessorImpl_GPU::InitKernel_max_pool_diff(const char* a_filePath)
{
  std::string shaderPath = AlterShaderPath("shaders_gpu/kernel3D_max_pool_diff.comp.spv"); 
  const VkSpecializationInfo* kspec = nullptr;
  max_pool_diffDSLayout = Createmax_pool_diffDSLayout();
  MakeComputePipelineAndLayout(shaderPath.c_str(), "main", kspec, max_pool_diffDSLayout, &max_pool_diffLayout, &max_pool_diffPipeline);
}


void TensorProcessorImpl_GPU::InitKernels(const char* a_filePath)
{
  InitKernel_cos(a_filePath);
  InitKernel_sin(a_filePath);
  InitKernel_urand(a_filePath);
  InitKernel_set_input(a_filePath);
  InitKernel_equal(a_filePath);
  InitKernel_min(a_filePath);
  InitKernel_log(a_filePath);
  InitKernel_logical_or(a_filePath);
  InitKernel_not_equal(a_filePath);
  InitKernel_sqrt(a_filePath);
  InitKernel_max(a_filePath);
  InitKernel_pow(a_filePath);
  InitKernel_greater(a_filePath);
  InitKernel_dilate(a_filePath);
  InitKernel_outer_product(a_filePath);
  InitKernel_pad(a_filePath);
  InitKernel_less(a_filePath);
  InitKernel_max_pool_3D_diff(a_filePath);
  InitKernel_less_equal(a_filePath);
  InitKernel_sum(a_filePath);
  InitKernel_copy(a_filePath);
  InitKernel_osum(a_filePath);
  InitKernel_not(a_filePath);
  InitKernel_exp(a_filePath);
  InitKernel_fill(a_filePath);
  InitKernel_get_output(a_filePath);
  InitKernel_logical_and(a_filePath);
  InitKernel_flip(a_filePath);
  InitKernel_sub(a_filePath);
  InitKernel_add(a_filePath);
  InitKernel_mul(a_filePath);
  InitKernel_greater_equal(a_filePath);
  InitKernel_div(a_filePath);
  InitKernel_minimum(a_filePath);
  InitKernel_max_pool_3D(a_filePath);
  InitKernel_maximum(a_filePath);
  InitKernel_transpose(a_filePath);
  InitKernel_matmul_transposed(a_filePath);
  InitKernel_max_pool(a_filePath);
  InitKernel_where(a_filePath);
  InitKernel_smax_diff(a_filePath);
  InitKernel_conv3d(a_filePath);
  InitKernel_conv2d(a_filePath);
  InitKernel_max_pool_diff(a_filePath);
}

void TensorProcessorImpl_GPU::InitBuffers(size_t a_maxThreadsCount, bool a_tempBuffersOverlay)
{
  ReserveEmptyVectors();
  
  m_maxThreadCount = a_maxThreadsCount;
  std::vector<VkBuffer> allBuffers;
  allBuffers.reserve(64);

  struct BufferReqPair
  {
    BufferReqPair() {  }
    BufferReqPair(VkBuffer a_buff, VkDevice a_dev) : buf(a_buff) { vkGetBufferMemoryRequirements(a_dev, a_buff, &req); }
    VkBuffer             buf = VK_NULL_HANDLE;
    VkMemoryRequirements req = {};
  };

  struct LocalBuffers
  {
    std::vector<BufferReqPair> bufs;
    size_t                     size = 0;
    std::vector<VkBuffer>      bufsClean;
  };

  std::vector<LocalBuffers> groups;
  groups.reserve(16);


  size_t largestIndex = 0;
  size_t largestSize  = 0;
  for(size_t i=0;i<groups.size();i++)
  {
    if(groups[i].size > largestSize)
    {
      largestIndex = i;
      largestSize  = groups[i].size;
    }
    groups[i].bufsClean.resize(groups[i].bufs.size());
    for(size_t j=0;j<groups[i].bufsClean.size();j++)
      groups[i].bufsClean[j] = groups[i].bufs[j].buf;
  }
  
  auto& allBuffersRef = allBuffers;

  m_classDataBuffer = vk_utils::createBuffer(device, sizeof(m_uboData),  VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | GetAdditionalFlagsForUBO());
  allBuffersRef.push_back(m_classDataBuffer);


  auto internalBuffersMem = AllocAndBind(allBuffersRef);
  if(a_tempBuffersOverlay)
  {
    for(size_t i=0;i<groups.size();i++)
      if(i != largestIndex)
        AssignBuffersToMemory(groups[i].bufsClean, internalBuffersMem.memObject);
  }
}

void TensorProcessorImpl_GPU::ReserveEmptyVectors()
{
  if(memory.capacity() == 0)
    memory.reserve(4);
}

void TensorProcessorImpl_GPU::InitMemberBuffers()
{
  std::vector<VkBuffer> memberVectors;
  std::vector<VkImage>  memberTextures;

  m_vdata.memoryBuffer = vk_utils::createBuffer(device, memory.capacity()*sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);
  memberVectors.push_back(m_vdata.memoryBuffer);


  AllocMemoryForMemberBuffersAndImages(memberVectors, memberTextures);
}




void TensorProcessorImpl_GPU::AssignBuffersToMemory(const std::vector<VkBuffer>& a_buffers, VkDeviceMemory a_mem)
{
  if(a_buffers.size() == 0 || a_mem == VK_NULL_HANDLE)
    return;

  std::vector<VkMemoryRequirements> memInfos(a_buffers.size());
  for(size_t i=0;i<memInfos.size();i++)
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkGetBufferMemoryRequirements(device, a_buffers[i], &memInfos[i]);
    else
    {
      memInfos[i] = memInfos[0];
      memInfos[i].size = 0;
    }
  }
  
  for(size_t i=1;i<memInfos.size();i++)
  {
    if(memInfos[i].memoryTypeBits != memInfos[0].memoryTypeBits)
    {
      std::cout << "[TensorProcessorImpl_GPU::AssignBuffersToMemory]: error, input buffers has different 'memReq.memoryTypeBits'" << std::endl;
      return;
    }
  }

  auto offsets = vk_utils::calculateMemOffsets(memInfos);
  for (size_t i = 0; i < memInfos.size(); i++)
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkBindBufferMemory(device, a_buffers[i], a_mem, offsets[i]);
  }
}

TensorProcessorImpl_GPU::MemLoc TensorProcessorImpl_GPU::AllocAndBind(const std::vector<VkBuffer>& a_buffers)
{
  MemLoc currLoc;
  if(a_buffers.size() > 0)
  {
    currLoc.memObject = vk_utils::allocateAndBindWithPadding(device, physicalDevice, a_buffers);
    currLoc.allocId   = m_allMems.size();
    m_allMems.push_back(currLoc);
  }
  return currLoc;
}

TensorProcessorImpl_GPU::MemLoc TensorProcessorImpl_GPU::AllocAndBind(const std::vector<VkImage>& a_images)
{
  MemLoc currLoc;
  if(a_images.size() > 0)
  {
    std::vector<VkMemoryRequirements> reqs(a_images.size()); 
    for(size_t i=0; i<reqs.size(); i++)
      vkGetImageMemoryRequirements(device, a_images[i], &reqs[i]);

    for(size_t i=0; i<reqs.size(); i++)
    {
      if(reqs[i].memoryTypeBits != reqs[0].memoryTypeBits)
      {
        std::cout << "TensorProcessorImpl_GPU::AllocAndBind(textures): memoryTypeBits warning, need to split mem allocation (override me)" << std::endl;
        break;
      }
    } 

    auto offsets  = vk_utils::calculateMemOffsets(reqs);
    auto memTotal = offsets[offsets.size() - 1];

    VkMemoryAllocateInfo allocateInfo = {};
    allocateInfo.sType           = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocateInfo.pNext           = nullptr;
    allocateInfo.allocationSize  = memTotal;
    allocateInfo.memoryTypeIndex = vk_utils::findMemoryType(reqs[0].memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, physicalDevice);
    VK_CHECK_RESULT(vkAllocateMemory(device, &allocateInfo, NULL, &currLoc.memObject));
    
    for(size_t i=0;i<a_images.size();i++) {
      VK_CHECK_RESULT(vkBindImageMemory(device, a_images[i], currLoc.memObject, offsets[i]));
    }

    currLoc.allocId = m_allMems.size();
    m_allMems.push_back(currLoc);
  }
  return currLoc;
}

void TensorProcessorImpl_GPU::FreeAllAllocations(std::vector<MemLoc>& a_memLoc)
{
  // in general you may check 'mem.allocId' for unique to be sure you dont free mem twice
  // for default implementation this is not needed
  for(auto mem : a_memLoc)
    vkFreeMemory(device, mem.memObject, nullptr);
  a_memLoc.resize(0);
}     

void TensorProcessorImpl_GPU::AllocMemoryForMemberBuffersAndImages(const std::vector<VkBuffer>& a_buffers, const std::vector<VkImage>& a_images)
{
  std::vector<VkMemoryRequirements> bufMemReqs(a_buffers.size()); // we must check that all buffers have same memoryTypeBits;
  for(size_t i = 0; i < a_buffers.size(); ++i)                    // if not, split to multiple allocations
  {
    if(a_buffers[i] != VK_NULL_HANDLE)
      vkGetBufferMemoryRequirements(device, a_buffers[i], &bufMemReqs[i]);
    else
    {
      bufMemReqs[i] = bufMemReqs[0];
      bufMemReqs[i].size = 0;
    }
  }

  bool needSplit = false;
  for(size_t i = 1; i < bufMemReqs.size(); ++i)
  {
    if(bufMemReqs[i].memoryTypeBits != bufMemReqs[0].memoryTypeBits)
    {
      needSplit = true;
      break;
    }
  }

  if(needSplit)
  {
    std::unordered_map<uint32_t, std::vector<uint32_t> > bufferSets;
    for(uint32_t j = 0; j < uint32_t(bufMemReqs.size()); ++j)
    {
      uint32_t key = uint32_t(bufMemReqs[j].memoryTypeBits);
      bufferSets[key].push_back(j);
    }

    for(const auto& buffGroup : bufferSets)
    {
      std::vector<VkBuffer> currGroup;
      for(auto id : buffGroup.second)
        currGroup.push_back(a_buffers[id]);
      AllocAndBind(currGroup);
    }
  }
  else
    AllocAndBind(a_buffers);

}


VkPhysicalDeviceFeatures2 TensorProcessorImpl_GPU::ListRequiredDeviceFeatures(std::vector<const char*>& deviceExtensions)
{
  static VkPhysicalDeviceFeatures2 features2 = {};
  features2.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
  features2.pNext = nullptr; 
  features2.features.shaderInt64   = false;
  features2.features.shaderFloat64 = false;  
  features2.features.shaderInt16   = false;
  void** ppNext = &features2.pNext;
  return features2;
}

TensorProcessorImpl_GPU::MegaKernelIsEnabled TensorProcessorImpl_GPU::m_megaKernelFlags;

