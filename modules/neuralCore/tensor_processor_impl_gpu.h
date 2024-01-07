#pragma once

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <array>

#include "vk_pipeline.h"
#include "vk_buffers.h"
#include "vk_utils.h"
#include "vk_copy.h"
#include "vk_context.h"


#include "tensor_processor_impl.h"

#include "include/TensorProcessorImpl_gpu_ubo.h"
class TensorProcessorImpl_GPU : public TensorProcessorImpl
{
public:

  TensorProcessorImpl_GPU() 
  {
  }
  virtual void InitVulkanObjects(VkDevice a_device, VkPhysicalDevice a_physicalDevice, size_t a_maxThreadsCount);
  
  virtual void SetVulkanContext(vk_utils::VulkanContext a_ctx) { m_ctx = a_ctx; }
  virtual void SetVulkanInOutFor_get_output(
    VkBuffer outBuffer,
    size_t   outOffset,
    uint32_t dummyArgument = 0)
  {
    get_output_local.outBuffer = outBuffer;
    get_output_local.outOffset = outOffset;
    InitAllGeneratedDescriptorSets_get_output();
  }

  virtual void SetVulkanInOutFor_set_input(
    VkBuffer inBuffer,
    size_t   inOffset,
    uint32_t dummyArgument = 0)
  {
    set_input_local.inBuffer = inBuffer;
    set_input_local.inOffset = inOffset;
    InitAllGeneratedDescriptorSets_set_input();
  }

  virtual void SetVulkanInOutFor_process(
    uint32_t dummyArgument = 0)
  {
    InitAllGeneratedDescriptorSets_process();
  }

  virtual ~TensorProcessorImpl_GPU();

  
  virtual void InitMemberBuffers();
  virtual void UpdateAll(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
  {
    UpdatePlainMembers(a_pCopyEngine);
    UpdateVectorMembers(a_pCopyEngine);
    UpdateTextureMembers(a_pCopyEngine);
  }
  
  virtual void CommitDeviceData(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyHelper) // you have to define this virtual function in the original imput class
  {
    InitMemberBuffers();
    UpdateAll(a_pCopyHelper);
  }  
  void CommitDeviceData() override { CommitDeviceData(m_ctx.pCopyHelper); }  
  void GetExecutionTime(const char* a_funcName, float a_out[4]) override; 
  

  virtual void ReserveEmptyVectors();
  virtual void UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  virtual void ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine);
  static VkPhysicalDeviceFeatures2 ListRequiredDeviceFeatures(std::vector<const char*>& deviceExtensions);
  
  virtual void get_outputCmd(VkCommandBuffer a_commandBuffer, float* out, unsigned offset, unsigned size);
  virtual void set_inputCmd(VkCommandBuffer a_commandBuffer, const float* in, unsigned offset, unsigned size);
  virtual void processCmd(VkCommandBuffer a_commandBuffer, const nn::TensorProgram &program);

  void get_output(float* out, unsigned offset, unsigned size) override;
  void set_input(const float* in, unsigned offset, unsigned size) override;
  void process(const nn::TensorProgram &program) override;

  inline vk_utils::ExecTime Getget_outputExecutionTime() const { return m_exTimeget_output; }
  inline vk_utils::ExecTime Getset_inputExecutionTime() const { return m_exTimeset_input; }
  inline vk_utils::ExecTime GetprocessExecutionTime() const { return m_exTimeprocess; }

  vk_utils::ExecTime m_exTimeget_output;
  vk_utils::ExecTime m_exTimeset_input;
  vk_utils::ExecTime m_exTimeprocess;

  virtual void copyKernelFloatCmd(uint32_t length);
  
  virtual void powCmd(float *data, unsigned steps, Variable A, Variable B, Variable C);
  virtual void expCmd(float *data, unsigned steps, Variable A, Variable B);
  virtual void sinCmd(float *data, unsigned steps, Variable A, Variable B);
  virtual void osumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);
  virtual void divCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C);
  virtual void transposeCmd(float *data, unsigned steps, unsigned row_len, unsigned col_len, Variable A, Variable B);
  virtual void addCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C);
  virtual void fillCmd(float *data, unsigned steps, Variable A, float val);
  virtual void logCmd(float *data, unsigned steps, Variable A, Variable B);
  virtual void subCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C);
  virtual void matmul_transposedCmd(float *data, unsigned A_row_len, unsigned A_col_len, unsigned B_col_len, 
                                          Variable A, Variable B, Variable C);
  virtual void outer_p_addCmd(float *data, unsigned step, unsigned A_len, unsigned B_len, 
                                         Variable A, Variable B, Variable C);
  virtual void get_outputCmd(float* data_out, unsigned offset, unsigned size);
  virtual void copyCmd(float *data, unsigned steps, unsigned from, unsigned to, Variable A, Variable B);
  virtual void set_inputCmd(const float* data_in, unsigned offset, unsigned size);
  virtual void sumCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);
  virtual void cosCmd(float *data, unsigned steps, Variable A, Variable B);
  virtual void mulCmd(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C);
  virtual void outer_productCmd(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
                                      Variable A, Variable B, Variable C);
  
  struct MemLoc
  {
    VkDeviceMemory memObject = VK_NULL_HANDLE;
    size_t         memOffset = 0;
    size_t         allocId   = 0;
  };

  virtual MemLoc AllocAndBind(const std::vector<VkBuffer>& a_buffers); ///< replace this function to apply custom allocator
  virtual MemLoc AllocAndBind(const std::vector<VkImage>& a_image);    ///< replace this function to apply custom allocator
  virtual void   FreeAllAllocations(std::vector<MemLoc>& a_memLoc);    ///< replace this function to apply custom allocator

protected:

  VkPhysicalDevice           physicalDevice = VK_NULL_HANDLE;
  VkDevice                   device         = VK_NULL_HANDLE;
  vk_utils::VulkanContext    m_ctx          = {};
  VkCommandBuffer            m_currCmdBuffer   = VK_NULL_HANDLE;
  uint32_t                   m_currThreadFlags = 0;
  std::vector<MemLoc>        m_allMems;
  VkPhysicalDeviceProperties m_devProps;

  VkBufferMemoryBarrier BarrierForClearFlags(VkBuffer a_buffer);
  VkBufferMemoryBarrier BarrierForSingleBuffer(VkBuffer a_buffer);
  void BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum);

  virtual void InitHelpers();
  virtual void InitBuffers(size_t a_maxThreadsCount, bool a_tempBuffersOverlay = true);
  virtual void InitKernels(const char* a_filePath);
  virtual void AllocateAllDescriptorSets();

  virtual void InitAllGeneratedDescriptorSets_get_output();
  virtual void InitAllGeneratedDescriptorSets_set_input();
  virtual void InitAllGeneratedDescriptorSets_process();

  virtual void AssignBuffersToMemory(const std::vector<VkBuffer>& a_buffers, VkDeviceMemory a_mem);

  virtual void AllocMemoryForMemberBuffersAndImages(const std::vector<VkBuffer>& a_buffers, const std::vector<VkImage>& a_image);
  virtual std::string AlterShaderPath(const char* in_shaderPath) { return std::string("") + std::string(in_shaderPath); }

  
  

  struct get_output_Data
  {
    VkBuffer outBuffer = VK_NULL_HANDLE;
    size_t   outOffset = 0;
    bool needToClearOutput = false;
  } get_output_local;

  struct set_input_Data
  {
    VkBuffer inBuffer = VK_NULL_HANDLE;
    size_t   inOffset = 0;
    bool needToClearOutput = false;
  } set_input_local;

  struct process_Data
  {
    bool needToClearOutput = false;
  } process_local;



  struct MembersDataGPU
  {
    VkBuffer memoryBuffer = VK_NULL_HANDLE;
    size_t   memoryOffset = 0;
  } m_vdata;
  
  
  size_t m_maxThreadCount = 0;
  VkBuffer m_classDataBuffer = VK_NULL_HANDLE;

  VkPipelineLayout      powLayout   = VK_NULL_HANDLE;
  VkPipeline            powPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout powDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatepowDSLayout();
  virtual void InitKernel_pow(const char* a_filePath);
  VkPipelineLayout      expLayout   = VK_NULL_HANDLE;
  VkPipeline            expPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout expDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateexpDSLayout();
  virtual void InitKernel_exp(const char* a_filePath);
  VkPipelineLayout      sinLayout   = VK_NULL_HANDLE;
  VkPipeline            sinPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout sinDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatesinDSLayout();
  virtual void InitKernel_sin(const char* a_filePath);
  VkPipelineLayout      osumLayout   = VK_NULL_HANDLE;
  VkPipeline            osumPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout osumDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateosumDSLayout();
  virtual void InitKernel_osum(const char* a_filePath);
  VkPipelineLayout      divLayout   = VK_NULL_HANDLE;
  VkPipeline            divPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout divDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatedivDSLayout();
  virtual void InitKernel_div(const char* a_filePath);
  VkPipelineLayout      transposeLayout   = VK_NULL_HANDLE;
  VkPipeline            transposePipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout transposeDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatetransposeDSLayout();
  virtual void InitKernel_transpose(const char* a_filePath);
  VkPipelineLayout      addLayout   = VK_NULL_HANDLE;
  VkPipeline            addPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout addDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreateaddDSLayout();
  virtual void InitKernel_add(const char* a_filePath);
  VkPipelineLayout      fillLayout   = VK_NULL_HANDLE;
  VkPipeline            fillPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout fillDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatefillDSLayout();
  virtual void InitKernel_fill(const char* a_filePath);
  VkPipelineLayout      logLayout   = VK_NULL_HANDLE;
  VkPipeline            logPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout logDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatelogDSLayout();
  virtual void InitKernel_log(const char* a_filePath);
  VkPipelineLayout      subLayout   = VK_NULL_HANDLE;
  VkPipeline            subPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout subDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatesubDSLayout();
  virtual void InitKernel_sub(const char* a_filePath);
  VkPipelineLayout      matmul_transposedLayout   = VK_NULL_HANDLE;
  VkPipeline            matmul_transposedPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout matmul_transposedDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout Creatematmul_transposedDSLayout();
  virtual void InitKernel_matmul_transposed(const char* a_filePath);
  VkPipelineLayout      outer_p_addLayout   = VK_NULL_HANDLE;
  VkPipeline            outer_p_addPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout outer_p_addDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout Createouter_p_addDSLayout();
  virtual void InitKernel_outer_p_add(const char* a_filePath);
  VkPipelineLayout      get_outputLayout   = VK_NULL_HANDLE;
  VkPipeline            get_outputPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout get_outputDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout Createget_outputDSLayout();
  virtual void InitKernel_get_output(const char* a_filePath);
  VkPipelineLayout      copyLayout   = VK_NULL_HANDLE;
  VkPipeline            copyPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout copyDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatecopyDSLayout();
  virtual void InitKernel_copy(const char* a_filePath);
  VkPipelineLayout      set_inputLayout   = VK_NULL_HANDLE;
  VkPipeline            set_inputPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout set_inputDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout Createset_inputDSLayout();
  virtual void InitKernel_set_input(const char* a_filePath);
  VkPipelineLayout      sumLayout   = VK_NULL_HANDLE;
  VkPipeline            sumPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout sumDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatesumDSLayout();
  virtual void InitKernel_sum(const char* a_filePath);
  VkPipelineLayout      cosLayout   = VK_NULL_HANDLE;
  VkPipeline            cosPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout cosDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatecosDSLayout();
  virtual void InitKernel_cos(const char* a_filePath);
  VkPipelineLayout      mulLayout   = VK_NULL_HANDLE;
  VkPipeline            mulPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout mulDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatemulDSLayout();
  virtual void InitKernel_mul(const char* a_filePath);
  VkPipelineLayout      outer_productLayout   = VK_NULL_HANDLE;
  VkPipeline            outer_productPipeline = VK_NULL_HANDLE; 
  VkDescriptorSetLayout outer_productDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout Createouter_productDSLayout();
  virtual void InitKernel_outer_product(const char* a_filePath);


  virtual VkBufferUsageFlags GetAdditionalFlagsForUBO() const;
  virtual uint32_t           GetDefaultMaxTextures() const;

  VkPipelineLayout      copyKernelFloatLayout   = VK_NULL_HANDLE;
  VkPipeline            copyKernelFloatPipeline = VK_NULL_HANDLE;
  VkDescriptorSetLayout copyKernelFloatDSLayout = VK_NULL_HANDLE;
  VkDescriptorSetLayout CreatecopyKernelFloatDSLayout();

  VkDescriptorPool m_dsPool = VK_NULL_HANDLE;
  VkDescriptorSet  m_allGeneratedDS[14];

  TensorProcessorImpl_GPU_UBO_Data m_uboData;
  
  constexpr static uint32_t MEMCPY_BLOCK_SIZE = 256;
  constexpr static uint32_t REDUCTION_BLOCK_SIZE = 256;

  virtual void MakeComputePipelineAndLayout(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, 
                                            VkPipelineLayout* pPipelineLayout, VkPipeline* pPipeline);
  virtual void MakeComputePipelineOnly(const char* a_shaderPath, const char* a_mainName, const VkSpecializationInfo *a_specInfo, const VkDescriptorSetLayout a_dsLayout, VkPipelineLayout pipelineLayout, 
                                       VkPipeline* pPipeline);

  std::vector<VkPipelineLayout> m_allCreatedPipelineLayouts; ///<! remenber them here to delete later
  std::vector<VkPipeline>       m_allCreatedPipelines;       ///<! remenber them here to delete later
public:

  struct MegaKernelIsEnabled
  {
    bool dummy = 0;
  };

  static MegaKernelIsEnabled  m_megaKernelFlags;
  static MegaKernelIsEnabled& EnabledPipelines() { return m_megaKernelFlags; }

};


