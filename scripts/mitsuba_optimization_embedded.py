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
                'filename': diffuse_tex_path
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

def init(base_path, image_w, image_h, spp, mitsuba_variant, render_style, texture_name):
  mi.set_variant(mitsuba_variant)
  scene_dict = {'type': 'scene'}
  scene_dict['integrator'] = {
    'type': 'direct_reparam',
    'reparam_rays': 8,
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
              'rfilter': { 'type': 'gaussian' },
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
    elif (render_style == "textured_const"):
      bsdf = {
          'type': 'diffuse',
                  'reflectance': {
                      "type": "bitmap",
                      "filename": base_path + "../textures/" + texture_name,
                      "filter_type": 'bilinear',
                      "wrap_mode": 'clamp'
                  },
      }
    else:
      print("Unknown render_style = ", render_style)

    scene_dict['model'] = {
        'type': 'obj',
        'filename': base_path + 'meshes/sphere.obj',
        'to_world': T.scale(0.67),
        'bsdf': bsdf
    }
    scene_dict['light'] = {
            'type': 'obj',
            'filename': base_path + 'meshes/sphere.obj',
            'emitter': {
                'type': 'area',
                'radiance': {'type': 'rgb', 'value': [100, 100, 100]}
            },
            'to_world': T.translate([0, 0.5, 10])
        }
    scene_dict['rectangle'] = {
          'type': 'obj',
          'filename': base_path + 'meshes/rectangle.obj',
          'to_world': T.translate([0, 0, -10]).scale(25),
          'emitter': {
              'type': 'area',
              'radiance': {'type': 'rgb', 'value': [1, 1, 1]}
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
    'spp' : spp
  }
  return context

def init_optimization(context, img_ref_dir, loss, save_intermediate_images):
  context['img_ref_dir'] = img_ref_dir
  context['loss_function'] = loss
  context['save_intermediate_images'] = int(save_intermediate_images)
  context['status'] = 'optimization_no_tex'

def init_optimization_with_tex(context, img_ref_dir, loss, save_intermediate_images):
  context['img_ref_dir'] = img_ref_dir
  context['loss_function'] = loss
  context['save_intermediate_images'] = int(save_intermediate_images)
  context['status'] = 'optimization_with_tex'

  opt = mi.ad.Adam(lr=0.1)
  opt['model.bsdf.reflectance.data'] = context['params']['model.bsdf.reflectance.data']
  context['tex_optimizer'] = opt

def render_and_save_to_file(context, save_filename):
  scene = context['scene']
  params = context['params']

  params['model.vertex_positions'] = context['vertex_positions']
  params['model.vertex_normals'] = context['vertex_normals']
  params['model.vertex_texcoords'] = context['vertex_texcoords']
  vertex_count = int(len(params['model.vertex_positions'])/3)
  params['model.vertex_count'] = vertex_count
  params['model.face_count'] = int(vertex_count/3)
  params['model.faces'] = list(range(vertex_count))
  params.update()

  img_ref = mi.render(scene, params, seed=0, spp=context['spp'])
  mi.util.write_bitmap(save_filename, img_ref)
  mi.util.convert_to_bitmap(img_ref)
  context['img_ref'] = img_ref
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

def render(it, context):
  if (int(it) == 0):
    with Image.open(context['img_ref_dir']) as img:
      img.load()
    img = img.convert('RGB')
    img_raw = numpy.asarray(img)
    img_raw = (img_raw/255.0) ** 2.2
    context['img_ref'] = mi.TensorXf(img_raw)
  scene = context['scene']
  params = context['params']
  img_ref = context['img_ref']

  params['model.vertex_positions'] = context['vertex_positions']
  params['model.vertex_normals'] = context['vertex_normals']
  params['model.vertex_texcoords'] = context['vertex_texcoords']
  vertex_count = int(len(params['model.vertex_positions'])/3)
  params['model.vertex_count'] = vertex_count
  params['model.face_count'] = int(vertex_count/3)
  if (len(params['model.faces']) != vertex_count):
    params['model.faces'] = list(range(vertex_count))
  params.update()

  dr.enable_grad(params['model.vertex_positions'])
  dr.enable_grad(params['model.vertex_normals'])
  dr.enable_grad(params['model.vertex_texcoords'])
  
  img = mi.render(scene, params, seed=it, spp=context['spp']) # image = F_render(scene)
  loss = context['loss_function'](img, img_ref) # loss = F_loss(image)
  dr.backward(loss)

  context['vertex_positions_grad'] = dr.grad(params['model.vertex_positions'])
  context['vertex_normals_grad'] = dr.grad(params['model.vertex_normals'])
  context['vertex_texcoords_grad'] = dr.grad(params['model.vertex_texcoords'])
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
      mi.util.write_bitmap("saves/tex_opt_iter"+str(it)+".png", mi.Bitmap(opt['model.bsdf.reflectance.data']))
  else:
    if (context['save_intermediate_images'] > 0 and int(it) % 10 == 0):
      mi.util.write_bitmap("saves/iter"+str(it)+".png", img)
      mi.util.write_bitmap("saves/iter"+str(it)+"_diff.png", dr.sqr(img - img_ref))
  return loss[0]

def get_params(context, key):
  lst = list(context[key])
  return struct.pack('%sf' % len(lst), *lst)

def set_params(context, key, bytes, n):
  tup = struct.unpack('%sf' % n, bytes)
  context[key] = tup