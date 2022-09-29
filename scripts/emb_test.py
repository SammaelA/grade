import drjit as dr
import mitsuba as mi
from matplotlib import pyplot as plt
import csv
from mitsuba.scalar_rgb import Transform4f as T
import struct

def init(base_path):
  image_size = 256
  print("base_path = ", base_path)
  mi.set_variant('cuda_ad_rgb')
  integrator = {
      'type': 'direct_reparam',
  }
  scene = mi.load_dict({
      'type': 'scene',
      'integrator': integrator,
      'sensor':  {
          'type': 'perspective',
          'to_world': T.look_at(
                          origin=(0, 0, 2),
                          target=(0, 0, 0),
                          up=(0, 1, 0)
                      ),
          'fov': 60,
          'film': {
              'type': 'hdrfilm',
              'width': image_size,
              'height': image_size,
              'rfilter': { 'type': 'gaussian' },
              'sample_border': True
          },
      },
      'model': {
          'type': 'obj',
          'filename': base_path + 'meshes/sphere.obj',
          'to_world': T.scale(0.67),
          'bsdf': {
              'type': 'diffuse',
              'reflectance': { 'type': 'rgb', 'value': (0.9, 0.8, 0.3) },
          },
      },
      'light': {
          'type': 'obj',
          'filename': base_path + 'meshes/sphere.obj',
          'emitter': {
              'type': 'area',
              'radiance': {'type': 'rgb', 'value': [1e3, 1e3, 1e3]}
          },
          'to_world': T.translate([2.5, 2.5, 7.0]).scale(0.25)
      }
  })

  #img_ref = mi.render(scene, seed=0, spp=256)
  #mi.util.write_bitmap("img_ref.png", img_ref)
  #mi.util.convert_to_bitmap(img_ref)
  params = mi.traverse(scene)
  context = {
    'scene' : scene,
    #'img_ref' : img_ref,
    'params' : params, 
    'vertex_positions' : params['model.vertex_positions'],
    'vertex_normals' : params['model.vertex_normals'],
    'vertex_texcoords' : params['model.vertex_texcoords'],
    'spp' : 16
  }
  return context

def get_sensor(image_w, image_h):
  sensor = {
          'type': 'perspective',
          'to_world': T.look_at(
                          origin=(0, 0, 2),
                          target=(0, 0, 0),
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
  return sensor

def init_optimization(context, image_w, image_h, spp, img_ref_dir):
  #context['params']['sensor'] = get_sensor(image_w, image_h)
  context['spp'] = spp
  #mi.util.write_bitmap

def render_and_save_to_file(context, image_w, image_h, spp, save_filename):
  #sensor_prev = context['params']['sensor']
  spp_prev = context['spp']

  #context['params']['sensor'] = get_sensor(image_w, image_h)
  context['spp'] = spp

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

  #context['params']['sensor'] = sensor_prev
  context['spp'] = spp_prev
def F_loss(img, img_ref):
    loss = dr.sqrt(dr.sum(dr.sqr(img - img_ref)) / len(img))
    return loss

def render(it, context):
  scene = context['scene']
  params = context['params']
  img_ref = context['img_ref']

  params['model.vertex_positions'] = context['vertex_positions']
  params['model.vertex_normals'] = context['vertex_normals']
  params['model.vertex_texcoords'] = context['vertex_texcoords']
  vertex_count = int(len(params['model.vertex_positions'])/3)
  params['model.vertex_count'] = vertex_count
  params['model.face_count'] = int(vertex_count/3)
  params['model.faces'] = list(range(vertex_count))
  params.update()

  dr.enable_grad(params['model.vertex_positions'])
  dr.enable_grad(params['model.vertex_normals'])
  dr.enable_grad(params['model.vertex_texcoords'])
  
  img = mi.render(scene, params, seed=it, spp=context['spp']) # image = F_render(scene)
  loss = F_loss(img, img_ref) # loss = F_loss(image)
  dr.backward(loss)

  context['vertex_positions_grad'] = dr.grad(params['model.vertex_positions'])
  context['vertex_normals_grad'] = dr.grad(params['model.vertex_normals'])
  context['vertex_texcoords_grad'] = dr.grad(params['model.vertex_texcoords'])

  if (int(it) % 100 == 0):
    mi.util.write_bitmap("saves/iter"+str(it)+".png", img)
  return loss[0]

def get_params(context, key):
  lst = list(context[key])
  return struct.pack('%sf' % len(lst), *lst)

def set_params(context, key, bytes, n):
  tup = struct.unpack('%sf' % n, bytes)
  context[key] = tup