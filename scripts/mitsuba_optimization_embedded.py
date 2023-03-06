import drjit as dr
import mitsuba as mi
from matplotlib import pyplot as plt
from PIL import Image
from mitsuba.scalar_rgb import Transform4f as T
import struct
import time
import numpy

def wood_material_principled(diffuse_tex_path: str, roughness_tex_path: str, specular=0.3):
    my_bsdf = mi.load_dict({
            'id': 'wood_material',
            'type': 'principled',
            'specular': specular,
            'base_color': {
                'type': 'bitmap',
                'filename': diffuse_tex_path
            },
            'roughness': {
                'type': 'bitmap',
                'filename': roughness_tex_path
            }
        })
    return my_bsdf


def metal_material_principled(diffuse_tex_path: str, roughness_tex_path: str, specular=0.85, metalness=0.85):
    my_bsdf = mi.load_dict({
            'id': 'metal_material',
            'type': 'principled',
            'specular': specular,
            'base_color': {
                'type': 'bitmap',
                'filename': diffuse_tex_path
            },
            'roughness': {
                'type': 'bitmap',
                'filename': roughness_tex_path
            },
            'metallic': metalness,
            'spec_tint': 1.0
        })
    return my_bsdf


def rough_conductor_cu(roughness_tex_path: str):
    my_bsdf = mi.load_dict({
            'id': 'metal_material',
            'type': 'roughconductor',
            'material': 'Cu',
            'distribution': 'ggx',
            'alpha':  {
                'type': 'bitmap',
                'filename': roughness_tex_path
            },
        })
    return my_bsdf


def porcelain_roughplastic(diffuse_tex_path: str, roughness=0.3):
    my_bsdf = mi.load_dict({
            'id': 'porcelain_material',
            'type': 'roughplastic',
            'distribution': 'ggx',
            'int_ior': 1.504,
            'diffuse_reflectance': {
                'type': 'bitmap',
                'filename': diffuse_tex_path,
                "filter_type": 'bilinear',
                "wrap_mode": 'clamp'
            },
            'alpha': roughness
            # 'alpha': {
            #     'type': 'bitmap',
            #     'filename': roughness_tex_path
            # }
        })
    return my_bsdf


def lambert(diffuse_tex_path: str):
    my_bsdf = mi.load_dict({
            'id': 'wood_material',
            'type': 'diffuse',
            'reflectance': {
                'type': 'bitmap',
                'filename': diffuse_tex_path
            }
        })
    return my_bsdf

def get_material_by_name(texture_name, material_name):
  if (material_name == "very smooth porcelain"):
    return porcelain_roughplastic(texture_name, 0)
  elif (material_name == "smooth porcelain"):
    return porcelain_roughplastic(texture_name, 0.001)
  elif (material_name == "porcelain"):
    return porcelain_roughplastic(texture_name, 0.01)
  elif (material_name == "ceramics"):
    return porcelain_roughplastic(texture_name, 0.1)
  elif (material_name == "rough ceramics"):
    return porcelain_roughplastic(texture_name, 0.3)
  else:
    print("unknown material name ", material_name)
    return porcelain_roughplastic(texture_name, 0)

def init(base_path, image_w, image_h, spp, mitsuba_variant, render_style, texture_name, material_name):
  texture_name = base_path + "../textures/" + texture_name
  mi.set_variant(mitsuba_variant)
  scene_dict = {'type': 'scene'}
  if (render_style == "monochrome_demo" or render_style == "textured_demo"):
    scene_dict['integrator'] = {
      'type': 'path',
      'max_depth': 8
    }
  elif (render_style == "silhouette"):
    scene_dict['integrator'] = {
      'type': 'emission_reparam',
      'reparam_rays': 8,
      'reparam_antithetic' : False
    }
  else:
    scene_dict['integrator'] = {
      'type': 'direct_reparam',
      'reparam_rays': 16,
      'reparam_antithetic' : False
    }
  scene_dict['sensor'] = {
          'type': 'perspective',
          'to_world': T.look_at(
                          origin=(0, 0.5, 1.5),
                          target=(0, 0.5, 0),
                          up=(0, 1, 0)
                      ),
          'fov': 60,
          'film': {
              'type': 'hdrfilm',
              'width': image_w,
              'height': image_h,
              'rfilter': { 'type': 'box'},
              'sample_border': True
          },
      }
  if (render_style == "silhouette"):
    scene_dict['model'] = {
          'type': 'obj',
          'filename': base_path + 'meshes/sphere.obj',
          'to_world': T.scale(0.67),
          'emitter': {
              'type': 'area',
              'radiance': {'type': 'rgb', 'value': [1, 1, 1]}
          }
    }
  else:
    bsdf = {
        'type': 'diffuse',
        'reflectance': {'type': 'rgb', 'value': (0.9, 0.8, 0.3)},
    }
    if (render_style == "monochrome"):
        bsdf = {
            'type': 'diffuse',
            'reflectance': {'type': 'rgb', 'value': (0.9, 0.8, 0.3)},
        }
    elif (render_style == "textured_const" or render_style == "textured_demo"):
      bsdf = get_material_by_name(texture_name, material_name)
    elif (render_style == "monochrome_demo"):
        bsdf = {
            'id': 'porcelain_material',
            'type': 'roughplastic',
            'distribution': 'ggx',
            'int_ior': 1.504,
            'diffuse_reflectance': {'type': 'rgb', 'value': (0.9, 0.8, 0.3)},
            'alpha': 0.001
        }
    else:
      print("Unknown render_style = ", render_style)

    scene_dict['model'] = {
        'type': 'obj',
        'filename': base_path + 'meshes/sphere.obj',
        'to_world': T.scale(0.67),
        "face_normals": True,
        'bsdf': bsdf
    }
    scene_dict['light'] = {
            'type': 'obj',
            'filename': base_path + 'meshes/sphere.obj',
            'emitter': {
                'type': 'area',
                'radiance': {'type': 'rgb', 'value': [100, 100, 100]}
            },
            'to_world': T.translate([0,0,0])
        }

    scene_dict['ambient_light'] = \
    {
      'type': 'constant',
      'radiance': {
          'type': 'rgb',
          'value': 0.2,
      }
    }


  scene = mi.load_dict(scene_dict)
  params = mi.traverse(scene)

  context = {
    'scene' : scene,
    'params' : params, 
    'vertex_positions' : params['model.vertex_positions'],
    'vertex_normals' : params['model.vertex_normals'],
    'vertex_texcoords' : params['model.vertex_texcoords'],
    'vertex_positions_grad' : params['model.vertex_positions'],
    'vertex_normals_grad' : params['model.vertex_normals'],
    'vertex_texcoords_grad' : params['model.vertex_texcoords'],
    'spp' : spp,
    'camera' : mi.load_dict(scene_dict['sensor']),
    'texture_name' : texture_name,
    'render_style' : render_style,
    'image_w' : image_w,
    'image_h' : image_h
  }
  if (render_style != "silhouette"):
    t = dr.unravel(mi.Point3f, params['light.vertex_positions'])
    context["light"] = dr.ravel(t)
  return context

def init_optimization(context, img_ref_dir, loss, learning_rate, cameras_count, save_intermediate_images):
  img_ref_dirs = img_ref_dir.split("#")
  for i in range(cameras_count):
    context['img_ref_dir_'+str(i)] = img_ref_dirs[i]
  context['loss_function'] = loss
  context['save_intermediate_images'] = int(save_intermediate_images)
  context['cameras_count'] = int(cameras_count)
  context['status'] = 'optimization_no_tex'

def init_optimization_with_tex(context, img_ref_dir, loss, learning_rate, cameras_count, save_intermediate_images):
  img_ref_dirs = img_ref_dir.split("#")
  for i in range(cameras_count):
    context['img_ref_dir_'+str(i)] = img_ref_dirs[i]
  context['loss_function'] = loss
  context['save_intermediate_images'] = int(save_intermediate_images)
  context['cameras_count'] = int(cameras_count)
  context['status'] = 'optimization_with_tex'

  opt = mi.ad.Adam(lr=learning_rate)
  opt['model.bsdf.diffuse_reflectance.data'] = context['params']['model.bsdf.diffuse_reflectance.data']
  context['tex_optimizer'] = opt

def set_camera(context, origin_x, origin_y, origin_z, target_x, target_y, target_z, up_x, up_y, up_z,
               fov, image_w, image_h):
  camera_dict = {
          'type': 'perspective',
          'to_world': T.look_at(
                          origin=(origin_x, origin_y, origin_z),
                          target=(target_x, target_y, target_z),
                          up=(up_x, up_y, up_z)
                      ),
          'fov': fov,
          'film': {
              'type': 'hdrfilm',
              'width': image_w,
              'height': image_h,
              'rfilter': { 'type': 'gaussian' },
              'sample_border': True
          },
    }
  context['camera'] =  mi.load_dict(camera_dict)  

def render_and_save_to_file(context, save_filename):
  scene = context['scene']
  params = context['params']

  params['model.vertex_positions'] = context['vertex_positions']
  vertex_count = int(len(params['model.vertex_positions'])/3)
  params['model.vertex_count'] = vertex_count
  params['model.face_count'] = int(vertex_count/3)
  params['model.faces'] = list(range(vertex_count))
  
  pos, angles, light_pos, ls_li_ali = get_scene_params(context)


  if (context["render_style"] != "silhouette"):
      light_transform = mi.Transform4f.translate([light_pos.x, light_pos.y, light_pos.z]).scale(ls_li_ali.x)
      t2 = dr.unravel(mi.Point3f, context["light"])
      lights_positions = light_transform @ t2
      params['light.vertex_positions'] = dr.ravel(lights_positions)
      params['light.emitter.radiance.value'] = mi.Color3f(ls_li_ali.y, ls_li_ali.y, ls_li_ali.y)

  if (context["render_style"] == "textured_const"):
      params['ambient_light.radiance.value'] = ls_li_ali.z

  t1 = dr.unravel(mi.Point3f, params['model.vertex_positions'])
  trafo = mi.Transform4f.translate([pos.x, pos.y, pos.z]).rotate([1, 0, 0], angles.x).rotate([0, 1, 0], angles.y).rotate([0, 0, 1], angles.z)
  tr_positions = trafo @ t1
  params['model.vertex_positions'] = dr.ravel(tr_positions)

  #we still should have proper size of normals and tc arrays
  if (len(params['model.vertex_positions']) != len(params['model.vertex_normals'])):
    params['model.vertex_normals'] = tuple([1] * len(params['model.vertex_positions']))
  if (2*len(params['model.vertex_positions']) != 3*len(params['model.vertex_texcoords'])):
    params['model.vertex_texcoords'] = tuple([1] * int(2*len(params['model.vertex_positions'])/3)) 
  
  params.update()
  t2_tensor = mi.Float(list(context['vertex_normals']))
  t2 = dr.unravel(mi.Normal3f, t2_tensor)
  trafo2 = mi.Transform4f.translate([pos.x, pos.y, pos.z]).rotate([1, 0, 0], angles.x).rotate([0, 1, 0], angles.y).rotate([0, 0, 1], angles.z)
  tr_normals = trafo2 @ t2
  context['vertex_normals'] = tuple(dr.ravel(tr_normals))
  params['model.vertex_normals'] = context['vertex_normals']
  params['model.vertex_texcoords'] = context['vertex_texcoords']

  max_samples = int(1024*1024*512 / (context['image_w']*context['image_h']))

  total_spp = context['spp']
  first = True
  while (total_spp > 0):
    spp = min(total_spp, max_samples)
    total_spp -= spp
    img = (float(spp)/context['spp']) * mi.render(scene, params, sensor = context['camera'], seed=0, spp=spp)
    if (first):
      img_ref = img
      first = False
    else:
      img_ref = img_ref + img
      
  mi.util.write_bitmap(save_filename, img_ref)

  time.sleep(1)

def F_loss_mse(img, img_ref):
    loss = dr.sum(dr.sqr(img - img_ref)) / len(img)
    return loss

def F_loss_mse_sqrt(img, img_ref):
    loss = dr.sqrt(dr.sum(dr.sqr(img - img_ref)) / len(img))
    return loss

def F_loss_mixed(img, img_ref):
    loss = dr.sum(0.5*dr.sqr(img - img_ref) + 0.5*dr.abs(img - img_ref)) / len(img)
    return loss

def get_scene_params(context):
    pos = mi.Point3f(context['camera_params'][0], 
                     context['camera_params'][1], 
                     context['camera_params'][2])
    tq = 180/numpy.pi
    angles = mi.Point3f(tq*context['camera_params'][3],
                        tq*context['camera_params'][4],
                        tq*context['camera_params'][5])
    
    light_pos = mi.Point3f(context['camera_params'][6], 
                           context['camera_params'][7], 
                           context['camera_params'][8])
    ls_li_ali   = mi.Point3f(context['camera_params'][9], #light_size, light_intensity, ambient_light_intensity
                           context['camera_params'][10], 
                           context['camera_params'][11])

    return (pos, angles, light_pos, ls_li_ali)
def render(it, context):

  #load reference on the first iteration
  if (int(it) == 0):
    for i in range(context['cameras_count']):
      print(context['img_ref_dir_'+str(i)])
      with Image.open(context['img_ref_dir_'+str(i)]) as img: ##TODO: proper link
        img.load()
      img = img.convert('RGB')
      img_raw = numpy.asarray(img)
      img_raw = (img_raw/255.0) ** 2.2
      if (context['status'] == 'optimization_with_tex'):
        img_flatten = img_raw.flatten()
        img_flatten_r = img_flatten[0::3]
        img_flatten_g = img_flatten[1::3]
        img_flatten_b = img_flatten[2::3]
        img_flatten_s = img_flatten_r + img_flatten_g + img_flatten_b
        img_flatten_s = numpy.repeat(img_flatten_s, 3)
        img_flatten[img_flatten_s > 1e-4] = 1
        img_flatten[img_flatten_s <= 1e-4] = 0
        img_mask = numpy.reshape(img_flatten, img_raw.shape)
      else:
        img_flatten = img_raw.flatten()
        img_mask = numpy.reshape(numpy.ones(img_flatten.shape, img_flatten.dtype), img_raw.shape)
      context['img_ref_mask_'+str(i)] = mi.TensorXf(img_mask)
      context['img_ref_'+str(i)] = mi.TensorXf(img_raw)
  
  #prepare model
  scene = context['scene']
  params = context['params']
  params['model.vertex_positions'] = context['vertex_positions']
  vertex_count = int(len(params['model.vertex_positions'])/3)
  params['model.vertex_count'] = vertex_count
  params['model.face_count'] = int(vertex_count/3)
  if (len(params['model.faces']) != vertex_count):
    params['model.faces'] = list(range(vertex_count))
  t1 = dr.unravel(mi.Point3f, params['model.vertex_positions'])
  camera_params_grad = []

  #render scene from each camera
  for camera_n in range(context['cameras_count']):
    img_ref = context['img_ref_'+str(camera_n)]
    img_ref_mask = context['img_ref_mask_'+str(camera_n)]

    pos, angles, light_pos, ls_li_ali = get_scene_params(context)

    dr.enable_grad(pos)
    dr.enable_grad(angles)
    dr.enable_grad(t1)

    if (context['status'] == 'optimization_with_tex'):
      dr.enable_grad(light_pos)
      dr.enable_grad(ls_li_ali)

      light_transform = mi.Transform4f.translate([light_pos.x, light_pos.y, light_pos.z]).scale(ls_li_ali.x)
      t2 = dr.unravel(mi.Point3f, context["light"])
      lights_positions = light_transform @ t2
      params['light.vertex_positions'] = dr.ravel(lights_positions)
      params['light.emitter.radiance.value'] = mi.Color3f(ls_li_ali.y, ls_li_ali.y, ls_li_ali.y)
      params['ambient_light.radiance.value'] = ls_li_ali.z

    trafo = mi.Transform4f.translate([pos.x, pos.y, pos.z]).rotate([1, 0, 0], angles.x).rotate([0, 1, 0], angles.y).rotate([0, 0, 1], angles.z)
    tr_positions = trafo @ t1
    params['model.vertex_positions'] = dr.ravel(tr_positions)
    
    #we still should have proper size of normals and tc arrays
    if (len(params['model.vertex_positions']) != len(params['model.vertex_normals'])):
      params['model.vertex_normals'] = tuple([1] * len(params['model.vertex_positions']))
    if (2*len(params['model.vertex_positions']) != 3*len(params['model.vertex_texcoords'])):
      params['model.vertex_texcoords'] = tuple([1] * int(2*len(params['model.vertex_positions'])/3)) 

    params.update()

    if (context['status'] == 'optimization_with_tex' or True):
      t2_tensor = mi.Float(list(context['vertex_normals']))
      t2 = dr.unravel(mi.Normal3f, t2_tensor)
      trafo2 = mi.Transform4f.translate([pos.x, pos.y, pos.z]).rotate([1, 0, 0], angles.x).rotate([0, 1, 0], angles.y).rotate([0, 0, 1], angles.z)
      tr_normals = trafo2 @ t2
      context['vertex_normals'] = tuple(dr.ravel(tr_normals))
      params['model.vertex_normals'] = context['vertex_normals']
      params['model.vertex_texcoords'] = context['vertex_texcoords']

    img = mi.render(scene, params, sensor = context['camera'], seed=it, spp=context['spp']) # image = F_render(scene)
    img = dr.minimum(dr.maximum(img, 0), 1)
    if (context['status'] == 'optimization_with_tex'):
      img = img_ref_mask * img
    loss = context['loss_function'](img, img_ref) # loss = F_loss(image)
    loss_PSNR = -10*(dr.log(1/loss)/dr.log(10)) #minus is because the pipeline tries to _minimize_ function
    dr.backward(loss_PSNR)

    camera_params_grad.append(dr.grad(pos).x)
    camera_params_grad.append(dr.grad(pos).y)
    camera_params_grad.append(dr.grad(pos).z)
    camera_params_grad.append(dr.grad(angles).x)
    camera_params_grad.append(dr.grad(angles).y)
    camera_params_grad.append(dr.grad(angles).z)
    if (context['status'] == 'optimization_with_tex'):
      camera_params_grad.append(dr.grad(light_pos).x)
      camera_params_grad.append(dr.grad(light_pos).y)
      camera_params_grad.append(dr.grad(light_pos).z)
      camera_params_grad.append(dr.grad(ls_li_ali).x)
      camera_params_grad.append(dr.grad(ls_li_ali).y)
      camera_params_grad.append(dr.grad(ls_li_ali).z)
    else:
      camera_params_grad.append([0])
      camera_params_grad.append([0])
      camera_params_grad.append([0])
      camera_params_grad.append([0])
      camera_params_grad.append([0])
      camera_params_grad.append([0])  

    if (camera_n == 0):
      vertex_positions_grad = dr.ravel(dr.grad(t1))
    else:
      vertex_positions_grad = vertex_positions_grad + dr.ravel(dr.grad(t1))
    
    if (context['save_intermediate_images'] > 0 and int(it) % 10 == 0):
      mi.util.write_bitmap("saves/iter"+str(it)+"_cam_"+str(camera_n)+".png", img)
      mi.util.write_bitmap("saves/iter"+str(it)+"_cam_"+str(camera_n)+"_diff.png", dr.sqr(img - img_ref))

  context['camera_params_grad'] = numpy.asarray(numpy.nan_to_num(camera_params_grad, nan=0, posinf=0, neginf=0))
  context['vertex_positions_grad'] = numpy.nan_to_num(vertex_positions_grad, nan=0, posinf=0, neginf=0) / float(context['cameras_count'])

  if (context['status'] == 'optimization_with_tex'):
    opt = context['tex_optimizer']
    opt.step()
    for key in opt.keys():
      if 'bsdf' in key:
        opt[key] = dr.clamp(opt[key], 0.0, 1.0)
      params[key] = opt[key]
    params.update()
    if context['save_intermediate_images'] > 0:
      mi.util.write_bitmap("saves/res_opt_iter"+str(it)+".png", img)
      mi.util.write_bitmap("saves/res_ref_opt_iter"+str(it)+".png", img_ref)
      mi.util.write_bitmap("saves/tex_opt_iter"+str(it)+".png", mi.Bitmap(opt['model.bsdf.diffuse_reflectance.data']))
    mi.util.write_bitmap(context['texture_name'], mi.Bitmap(opt['model.bsdf.diffuse_reflectance.data']))
  return loss_PSNR[0]

def get_params(context, key):
  lst = list(context[key])
  return struct.pack('%sf' % len(lst), *lst)

def set_params(context, key, bytes, n):
  tup = struct.unpack('%sf' % n, bytes)
  context[key] = tup