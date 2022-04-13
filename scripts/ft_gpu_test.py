import tensorflow as tf
tf.config.list_physical_devices('GPU')
sys_details = tf.sysconfig.get_build_info()
cuda_version = sys_details["cuda_version"]
print(cuda_version)
print(sys_details)
print(tf.sysconfig.get_include())